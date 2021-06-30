% Load and match image with the given sample jpg
clear all;clc;close all;
%Set the path to the desired data directory in MSTAR CD data
path2mstar='E:/BoxS/Box Sync/mstar/Data_from_CD/CD_2_publicTargets/TARGETS/TRAIN/17_DEG';
fileNamePrefix = {'2S1','BMP2_SN_9563','BMP2_SN_9566','BMP2_SN_C21',...
                'BRDM_2','BTR_60','BTR70_SN_C71','D7','T62',...
                'T72_SN_132','T72_SN_812','T72_SN_S7','ZIL131','ZSU_23_4'};
extensions = {'000','000','001','002','001','003','004','005','016',...
            '015','016','017','025','026'};
        
%Define constant parameters
f_center = 9.6e9;
velLight = 299792458;
numPixelsCrop = 128;
numPixelsTarget = 100;
bandwidth = 521e6;
fLower = f_center - bandwidth/2;
taylorWindow = kron(taylorwin(numPixelsTarget,4,-35),taylorwin(numPixelsTarget,4,-35).');
f = linspace(fLower,fLower + bandwidth,numPixelsTarget ).';
thetas = (-1.5:0.03:1.5-0.03);
fRep  = repmat(f,1,numPixelsTarget);

%% Iterate over all classes/folders

for idxPrefix = 1%:length(fileNamePrefix)
    fprintf('Processing Class %s ...\n',fileNamePrefix{idxPrefix});
    pathLoad=sprintf('%s/%s/',path2mstar,fileNamePrefix{idxPrefix});
    file_names=dir([pathLoad sprintf('*.%s',extensions{idxPrefix})]);
    file_names=cell2mat((extractfield(file_names,'name'))');
    
    %initialize the arrays to store
    lenth=length(file_names);
    arr_azi = zeros(1,lenth);
    depression=zeros(1,lenth);
    arr_img_comp=zeros(lenth,numPixelsCrop,numPixelsCrop);
    arr_img_fft_polar = zeros(numPixelsTarget,numPixelsTarget,lenth);

    %Iterate over all samples
    for i = 1:lenth
        path1=[pathLoad file_names(i,:)];
        gg=MSTAR_LOAD_IMAGE(path1);
        
        img_comp=(flipud(gg.ImageData));
        azi=gg.TargetAz;
        depression(i) = gg.MeasuredDepression;
        arr_azi(i)=azi;
        
        %"Square" the image
        if i==1
            len=min(size(img_comp));
        end
        img_comp=img_comp(1:len,1:len);
        
        %Crop to fixed size numPixelsCrop while maintaining centering
        centerIm = floor(size(img_comp)/2);
        img_comp = img_comp(centerIm(1) -floor(numPixelsCrop/2)+1:centerIm(1) +floor(numPixelsCrop/2),...
             centerIm(2) -floor(numPixelsCrop/2)+1:centerIm(2) +floor(numPixelsCrop/2));
        centerIm = floor(size(img_comp)/2);
        arr_img_comp(i,:,:) = img_comp;
        
        %form the polar meshgrid
        thetaRep = repmat(thetas,numPixelsTarget,1);
        k_1 =  (4*pi/velLight*cosd(depression(i))*fRep.*sind(thetaRep));
        k_2 =  (4*pi/velLight*cosd(depression(i))*fRep.*cosd(thetaRep));
        yy = linspace(min(k_2(:)),max(k_2(:)),100);
        xx = linspace(min(k_1(:)),max(k_1(:)),100);
        [XX,YY] = meshgrid(xx,yy); 
        
        %form image-domain cropping mask of desired size (64-by-64)
        mask1 = zeros(numPixelsCrop,numPixelsCrop);
        mask1(centerIm-32:centerIm+31,centerIm-32:centerIm+31) = 1;
        
        %Mask and transform to phase history domain
        arr_img_fft = fftshift(fft2(ifftshift(mask1.*img_comp)));
        
        %Undo the taylor window and crop to numPixelsTarget
        arr_img_fft_crop =(1./taylorWindow.*arr_img_fft(centerIm(1)-round(numPixelsTarget/2)+1:centerIm(1)...
            + numPixelsTarget -round(numPixelsTarget/2) ,centerIm(2)-round(numPixelsTarget/2)+1:centerIm(2) ...
        + numPixelsTarget -round(numPixelsTarget/2),1));
        arr_img_fft_polar(:,:,i) = (interp2(XX,YY,(arr_img_fft_crop),k_1,k_2,'spline',0));
    end
    
    %Save data to disc
    save(sprintf('../data/phase_histories/%s_PH',fileNamePrefix{idxPrefix}),...
          'arr_azi','arr_img_fft_polar','arr_img_comp','depression');
end