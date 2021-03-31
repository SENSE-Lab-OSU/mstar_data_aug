clc;
clear all;
close all;
%gpuDevice(1);d=gpuDevice;d.FreeMemory;

addpath('mtimesx')
path2PH=('../data/phase_histories/');
path2coeff=('../data/recovered_coefficients/');
sample_start=1; %Change sample_start to start generating from that sample
sample_end=1; %This is mostly dummy variable unless line 40 is commented
%Rough time calc: 2747*4*22=67 vs 2747*4*11=67
%Estimated time per sample (~80gb RAM)=  255.87-153.62= 102.25
%Estimated time per sample (~30gb RAM)=  237.87-136.3= 101.57
%figure;imagesc(abs(reshape(imgTrain(4,:,:),64,64)));

fileNamePrefix = {'BMP2_SN_9563', 'BTR_70_SN_C71', 'T72_SN_132', 'BTR_60','2S1', 'BRDM_2', 'D7', 'T62', 'ZIL131', 'ZSU_234'};
%load('red_idx_array');
shiftsy =[0,0.15,0, 0.15];
shiftsx =[0,0, 0.15,0.15];
array_dtheta=[-6:0.25:-0.25 0.25:0.25:6];
tol_add=[-ones(1,9)*0.01 -0.03 0 -0.01 -0.03 0 0 -0.02 zeros(1,3) 0.02 ...
    zeros(1,2) 0.02*ones(1,4) zeros(1,2) 0.02 zeros(1,6) -0.02 ...
    zeros(1,2) -0.02 0 -0.02 -0.03 zeros(1,2) -0.02 0 -0.02 -0.03];
tol_add=flip(-tol_add); %New change after dtheta sign fiasco!

for idxClass = 3%:length(fileNamePrefix)
    %e=load(sprintf('%s%s_trainImPhaseHistories',path2PH,fileNamePrefix{idxClass}));
    %e=load(sprintf('%s%s_masked',path2PH,fileNamePrefix{idxClass}));
    e=load(sprintf('%s%s_PH',path2PH,fileNamePrefix{idxClass}));

    
    m=load(sprintf('%s%s',path2coeff,fileNamePrefix{idxClass}));
    class_folder=sprintf('../data/gen_aug_data/%s',fileNamePrefix{idxClass});
    if ~isfolder(class_folder)
       mkdir(class_folder)
    end
    numRangeBins = 100;

    
    numTrainingSamples = size(m.y_recovered,3);
    disp('Line 38 will over-ride sample_end to iterate over all remaining samples');
    %sample_end=numTrainingSamples;
    
    taylorWindow = kron(taylorwin(100,4,-35),taylorwin(100,4,-35).');
    %thetas = [-5.5:0.03:-1.5-0.03 linspace(-1.5,1.5,100) 1.53:0.03:5.5 ] ;
    thetas = [-7.5:0.03:-1.5-0.03 linspace(-1.5,1.5,100) 1.53:0.03:7.5 ] ;

    bisectorAngles = 90 + thetas;
    azimuthVals = bisectorAngles;
    
    
    
    numFreqBins = 100;
    f_center = 9.6e9;
    bandwidth = 521e6;
    delF = bandwidth/100;
    fLower = f_center - bandwidth/2;
    f = linspace(fLower,fLower + bandwidth,100 ).';
    
    
    AzimuthBasisCenterSpacing = 0.4; % Degrees
    numAzimuthBasisCenters = round(3/AzimuthBasisCenterSpacing);
  
    
    
    
    %THe angular and spatial grid points
    
    L=30;
    azimuthBasisCenters = 90+ linspace(-1.5,1.5,numAzimuthBasisCenters);
    
    
    bisectorRep = repmat(bisectorAngles,length(f),1);
    
    
    xGrids = -L/2:0.3:L/2-0.3;
    yGrids = -L/2:0.3:L/2-0.3;
    
    [X,Y] = meshgrid(xGrids,yGrids);
    Xp = repmat(X(:)',numFreqBins,1);
    Yp = repmat(Y(:)',numFreqBins,1);
    X = X(:);
    Y = Y(:);
    F = repmat(f,1,numRangeBins^2);
    
    thetas1 = linspace(-1.5,1.5,100);
    fRep  = repmat(f,1,numFreqBins);
    
    velLight = 3e8;
    
    % Form Matrix A_mod
    % %Combined Implementation using MATLAB broadcasting reduced time but avg.
    
%     elevation=mean(e.depression);%Approximation
%     %stad=std(e.depression);
%     MC=1i*4*pi*F*cosd(elevation)/velLight; % 100,10000
%     MX=(Xp+ reshape(shiftsx,1,1,1,[])).*reshape(cosd(azimuthVals),1,1,[]); %100,10000,366,4
%     MY=(Yp+ reshape(shiftsy,1,1,1,[])).*reshape(sind(azimuthVals),1,1,[]); %100,10000,366,4
%     A_mod=1/sqrt(numFreqBins)*exp((MX+MY).*MC); %100,10000,366,4
    
    rangeResolution =  velLight/(2*bandwidth);
    BasisFuncTemp=[];
    bb = bisectorAngles(((bisectorAngles >=90-1.5 ) & (bisectorAngles <= 90+1.5)));
    
    numAzimuthBinsTotal = length(azimuthVals);
    pulseSelectIDx = 1:length(azimuthVals);
    

    for idxTrain = sample_start:sample_end
        gaussWidth = m.gaussWidthStore(idxTrain);
        dist = pdist2(bb.',azimuthBasisCenters');
        BasisFuncTemp = exp(-0.5*dist.^2/gaussWidth^2);
        normBasisFunc = diag((sum(BasisFuncTemp.^2,1)).^0.5);
        
        dist = pdist2(azimuthVals.',azimuthBasisCenters');
        BasisFunc = exp(-0.5*dist.^2/gaussWidth^2);
        BasisFunc = BasisFunc/normBasisFunc;
        
        %elevation = e.depression(idxTrain);
        elevation = e.elevation(idxTrain);

        
        %Using this discrete approach, we could possibly save full 128X128 
        %images
        imgTrain =zeros((length(array_dtheta)+1)*length(shiftsy),64,64);
        aziTrain=zeros((length(array_dtheta)+1)*length(shiftsy),1);
        elev=zeros((length(array_dtheta)+1)*length(shiftsy),1);
        countIm=1;
        
        for idxShifts =1:length(shiftsy)
            fprintf('processing class %d, image %d, shift %d \n'...
                ,idxClass,idxTrain,idxShifts);
            fNew = linspace(fLower*sind(90 + 1.5),fLower + bandwidth,100 ).';
            yy =   (4*pi/velLight*(fNew)*cosd(elevation));% (4*pi/velLight*f*cosd(elevation)); %
            xx1 = thetas1.';
            xx = (4*pi/velLight*f(end)*cosd(elevation)*cosd(90+xx1));
            [XX,YY] = meshgrid(xx,yy);
            pointsOrig = [XX(:).';YY(:).'];
            
            
            
            XX = reshape(XX,100,100);
            YY = reshape(YY,100,100);
            
            thetaRep = repmat(thetas1,100,1);
            
            
            k_1 =  (4*pi/velLight*cosd(elevation)*(fRep).*cosd(90+thetaRep));
            k_2 =  (4*pi/velLight*cosd(elevation)*(fRep).*sind(90+thetaRep));
            
            
            
            
            %Implementation using MATLAB broadcasting reduced time by 1/2  
            MC=1i*4*pi*F*cosd(elevation)/velLight;
            MX=(Xp+ shiftsx(idxShifts)).*reshape(cosd(azimuthVals),1,1,[]);
            MY=(Yp+ shiftsy(idxShifts)).*reshape(sind(azimuthVals),1,1,[]);
            A_mod=1/sqrt(numFreqBins)*exp((MX+MY).*MC);
            

%             % Added computation on GPU made it 1/3
%             d.FreeMemory
%             MC=gpuArray(1i*4*pi*F*cosd(elevation)/velLight);
%             MXY=(Xp+ shiftsx(idxShifts)).*reshape(cosd(azimuthVals),1,1,[]);
%             MXY=gpuArray(MXY+(Yp+ shiftsy(idxShifts)).*reshape(sind(azimuthVals),1,1,[]));
%             A_mod=1/sqrt(numFreqBins)*gather(exp(MXY.*MC));
            
            %Check for rounding errors
            %diff_mat=abs(A_mod-A_mod2);
            %aa=diff_mat(:,:,2);
            
            
            A = @(x,mode) SAR_operator_gen(x,mode,pulseSelectIDx,...
                numRangeBins,numFreqBins,azimuthVals,BasisFunc,...
                A_mod(:,:,pulseSelectIDx));

            fRep1= repmat(f,1,length(azimuthVals));
            y_recovered1 =   reshape(A(m.x_recovered(:,idxTrain),1),100,length(azimuthVals));
            
            % 0 degree case
            %y =  (y_recovered1(:,find((bisectorAngles >=90-1.5 ) & (bisectorAngles <= 90+1.5)))) ;
            y =  (y_recovered1(:,((bisectorAngles >=90-1.5) & (bisectorAngles <= 90+1.5)))) ;

            yTotal = y;
            ss = scatteredInterpolant(k_1(:),k_2(:),yTotal(:),'natural','nearest');
            arr_img_fft_cart1 = reshape(ss(XX(:),YY(:)),100,100);
            arr_img_fft_cart = zeros(128,128);
            arr_img_fft_cart(14:113,14:113) =   taylorWindow.*arr_img_fft_cart1 ;
            arr_img_recons = ifftshift(ifft2(arr_img_fft_cart));
            
            
            %yTotal= exp(1i*4*pi*fRep*cosd(elevation)/velLight.*...
            %    ( shiftsy(idxShifts).*sind(90+thetaRep) )).*e.arr_img_fft_polar(:,:,idxTrain);
            yTotal= exp(1i*4*pi*fRep*cosd(elevation)/velLight.*...
                ( shiftsy(idxShifts).*sind(90+thetaRep) +shiftsx(idxShifts).*cosd(90+thetaRep) ))...
                .*e.arr_img_fft_polar(:,:,idxTrain);
            normT= norm(yTotal,'fro');
            ss = scatteredInterpolant(k_1(:),k_2(:),yTotal(:),'natural','nearest');
            arr_img_fft_cart1 = reshape(ss(XX(:),YY(:)),100,100);
            arr_img_fft_cart = zeros(128,128);
            arr_img_fft_cart(14:113,14:113) =   taylorWindow.*arr_img_fft_cart1 ;
            arr_img_original = ifftshift(ifft2(arr_img_fft_cart));
            arr_img_residual = arr_img_original - arr_img_recons;
            imgTrain(countIm,:,:) = arr_img_original(32:95,32:95);
            aziTrain(countIm) = e.arr_azi(idxTrain);
            elev(countIm) = elevation;
            countIm = countIm+1;
            
            for idxtheta=1:length(array_dtheta)
                %the 2nd dtheta in next line seems to have some neigborhood
                %values at times
                dtheta=array_dtheta(idxtheta);
                %Very Important change!! Angle of rotation seemed reverse
                %so changed the sign only and reversed tol_add
                
                %y = (y_recovered1(:,((bisectorAngles >=90-1.5+dtheta) & ...
                %    (bisectorAngles <= 90+1.5+(dtheta+tol_add(1,idxtheta))))));
                y = (y_recovered1(:,((bisectorAngles >=90-1.5-dtheta) & ...
                    (bisectorAngles <= 90+1.5-(dtheta+tol_add(1,idxtheta))))));
                yTotal = normT*y/norm(y,'fro') ;
                ss = scatteredInterpolant(k_1(:),k_2(:),yTotal(:),'natural','nearest');
                arr_img_fft_cart1 = reshape(ss(XX(:),YY(:)),100,100);
                arr_img_fft_cart = zeros(128,128);
                arr_img_fft_cart(14:113,14:113) =   taylorWindow.*arr_img_fft_cart1 ;
                arr_img_recons = ifftshift(ifft2(arr_img_fft_cart));
                %arr_img_recons= imrotate(arr_img_recons,dtheta,'crop');
                %arr_img_recons = arr_img_recons +imrotate(arr_img_residual,-dtheta,'crop');
                arr_img_recons = arr_img_recons +imrotate(arr_img_residual,dtheta,'crop');

                imgTrain(countIm,:,:) = arr_img_recons(32:95,32:95);

                aziTrain(countIm) = e.arr_azi(idxTrain)+dtheta;
                elev(countIm) = elevation;
                countIm = countIm+1;
            end
            %save(sprintf('gen_aug_data\\%s_corrected_shifts_new',fileNamePrefix{idxClass}),'imgTrain','aziTrain','elev','-v7.3');
            %save(sprintf('%s\\sample_%d',class_folder,idxTrain),'imgTrain','aziTrain','elev');
        end
        %Only save once at the end of all subpixel shifts of the sample
        save(sprintf('%s\\sample_%d',class_folder,idxTrain),'imgTrain','aziTrain','elev');
    end
end