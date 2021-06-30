clc;
clear all;
close all;

addpath('spgl1-2.1');
addpath('mtimesx');

fileNamePrefix = {'2S1','BMP2_SN_9563','BMP2_SN_9566','BMP2_SN_C21',...
                'BRDM_2','BTR_60','BTR70_SN_C71','D7','T62',...
                'T72_SN_132','T72_SN_812','T72_SN_S7','ZIL131','ZSU_23_4'};
gaussWidth = 1;

%% Iterate over all classes

for idxClass =1:length(fileNamePrefix)
    data_PH=load(sprintf('../data/phase_histories/%s_PH',fileNamePrefix{idxClass}));
    numImages = size(data_PH.arr_img_fft_polar,3);
    f_center = 9.6e9;
    bandwidth = 521e6;
    delF = bandwidth/100;
    fLower = f_center - bandwidth/2;
    f = linspace(fLower,fLower + bandwidth,100 ).';
    
    % Using the horizontal polarization measurements
    thetas = (-1.5:0.03:1.5-0.03);
    bisectorAngles = 90 + thetas;
    azimuthVals = bisectorAngles;
    deltaAngles = 0;
    numChips = size(data_PH.arr_img_fft_polar,3);   
    velLight = 299792458;
    
    numFreqBins = length(f);
    pixelResolutionMSTAR = 0.202;
    numRangeBins = numFreqBins; % Number of range bins are calculated from the range resolution
    
    numAzimuthBinsTotal = length(bisectorAngles); % Number of azimuth looks are calculated across the range [0, 3) with resolution dependent on cross-range extent
    AzimuthBasisCenterSpacing = 0.4; % Degrees
    numAzimuthBasisCenters = round(3/AzimuthBasisCenterSpacing);
    
    gaussWidthMax=5;
    gaussWidthMin = 1;
    
    %The angular and spatial grid points
    L =30;
    azimuthBasisCenters = 90+ linspace(-1.5,1.5,numAzimuthBasisCenters);

    %using range resolution=0.3
    xGrids = -L/2:0.3:L/2-0.3;
    yGrids = -L/2:0.3:L/2-0.3;
    
    [X,Y] = meshgrid(xGrids,yGrids);
    Xp = repmat(X(:)',numFreqBins,1);
    Yp = repmat(Y(:)',numFreqBins,1);
    X = X(:);
    Y = Y(:);
    F = repmat(f,1,numRangeBins^2);
    
    BasisFuncTemp=[];
    groups_l12 = [];
    
    dist = pdist2(azimuthVals.',azimuthBasisCenters');
    BasisFunc = exp(-0.5*dist.^2/gaussWidth^2);
    normBasisFunc = diag((sum(BasisFunc.^2,1)).^0.5);
    BasisFunc = BasisFunc/normBasisFunc;
    numVariables= size(BasisFunc,2);
    groups_l12 = repmat((1:numRangeBins^2).',1,numVariables);
    groups_l12=groups_l12(:);
    
    numPulses = numAzimuthBinsTotal;
    pulseSelectIDx = 1:length(azimuthVals);
    az =azimuthVals;
    
    
    % loop over different chips
    y_recovered = zeros(100,100,numChips);
    y_residual =  zeros(100,100,numChips);
    x_recovered = zeros(100*100*numVariables,numChips);
    fileName = sprintf('../data/recovered_coefficients/%s',fileNamePrefix{idxClass});
    gaussWidthStore=zeros(numChips,1);
    for idxChips = 1%:numChips
        fprintf('processing class=%d,image=%d\n',idxClass,idxChips);
        
        depression = data_PH.depression(idxChips);
        A_mod = zeros(numFreqBins,numRangeBins^2,numAzimuthBinsTotal);
        for idxPulses=1:numAzimuthBinsTotal
            A_mod(:,:,idxPulses) =1/sqrt(numFreqBins)*exp(1i*4*pi*F*...
                cosd(depression)/velLight.*(Xp.*cosd(azimuthVals(idxPulses)) +...
                Yp.*sind(azimuthVals(idxPulses)) ));
        end
        
        
        y1= (data_PH.arr_img_fft_polar(:,:,idxChips));
        y1=y1(:); 
        snr = 20;
        numSamples = length(y1);
        sigma_n = sqrt(2)*norm(y1)*10^(-snr/20);
        
        options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',4);
        
        ff = @(g) reconError(y1,g,pulseSelectIDx,dist,azimuthVals,numFreqBins,numRangeBins,A_mod,groups_l12,sigma_n);
        tic;
        gaussWidth = fminbnd(ff,gaussWidthMin,gaussWidthMax,options);
        toc;
        
        %TODO: Looks redundant as in reconError.
        BasisFunc = exp(-0.5*dist.^2/gaussWidth^2);
        normBasisFunc = diag((sum(BasisFunc.^2,1)).^0.5);
        BasisFunc = BasisFunc/normBasisFunc;
        A = @(x,mode) SAR_operator_gen(x,mode,pulseSelectIDx,numRangeBins,numFreqBins,azimuthVals,BasisFunc,A_mod(:,:,pulseSelectIDx));
        opts = spgSetParms('iscomplex',1,'verbosity',0);
        %Find Optimum C
        C = spg_group(A, y1, groups_l12, sigma_n, opts );
        
        %Get Derived quantities
        x_recovered(:,idxChips) = C;
        y_recovered(:,:,idxChips) = reshape(A(C,1),100,100);
        y_residual(:,:,idxChips) = fliplr(data_PH.arr_img_fft_polar(:,:,idxChips)) - y_recovered(:,:,idxChips);
        gaussWidthStore(idxChips) = gaussWidth;
    end
    save(fileName,'x_recovered','y_recovered','y_residual','gaussWidthStore');
end