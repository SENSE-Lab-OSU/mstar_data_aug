% Load and match image with the given sample jpg
clear all;close all;clc;
% TODO: Specify path to the specific class' data from CD
pathLoad='E:/BoxS/Box Sync/mstar/Data_from_CD/CD_2_publicTargets/TARGETS/TRAIN/17_DEG/T72_SN_132/';
file_ext='015';
file_out_prefix = 'T72_SN_132';
% TODO: Specify path to the Solver_fixedBistatic_Group
%addpath('E:/BoxS/Box Sync/mstar/Solver_fixedBistatic_Group/');

% Load data
file_names=dir([pathLoad '*.' file_ext]);
file_names=cell2mat((extractfield(file_names,'name'))');

% Set some predefined parameters specific to processing MSTAR dataset
taylorWindow = kron(taylorwin(100,4,-35),taylorwin(100,4,-35).');
f_center = 9.6e9;
bandwidth = 521e6;
delF = bandwidth/100;
fLower = f_center - bandwidth/2;
f = linspace(fLower,fLower + bandwidth,100 ).';
velLight = 3e8;
fRep = repmat(f,1,100);
thetas = (-1.5:0.03:1.5-0.03);
numFreqs  = length(f);
numBisectors = length(thetas);
fSimulationRep = repmat(f,1,numBisectors);

%{
%Preallocate arrays
n_samples=size(file_names,1)
arr_azi=zeros(
arr_img_comp=
arr_img_comp_rot=
arr_img_fft=
arr_img_fft_crop=
arr_img_fft_polar=
%}

for i = 1:size(file_names,1)
    path1=[pathLoad file_names(i,:)];
    gg=MSTAR_LOAD_IMAGE(path1);
    elevation = gg.MeasuredDepression;
    yy = (4*pi/velLight*f*cosd(elevation));
    xx1 = (-1.5:0.03:1.5-0.03).';
    xx = (4*pi/velLight*f(end)*cosd(elevation)*sind(xx1));
    [XX,YY] = meshgrid(xx,yy);
    pointsOrig = [XX(:).';YY(:).'];
    img_comp=(flipud(gg.ImageData));
    azi=gg.TargetAz;
    rotationMat = [cosd(0*azi) -sind(0*azi);sind(0*azi) cosd(0*azi)];
    pointsRot=  rotationMat *pointsOrig;
    XX  = pointsRot(1,:);
    YY = pointsRot(2,:);
    XX = reshape(XX,100,100);
    YY = reshape(YY,100,100);
    
    len=128;
    img_comp=img_comp(1:len,1:len);
    img_comp_rot = imrotate(img_comp,azi,'bilinear','crop');

    arr_azi(i)=azi;
    
    thetasActual = 90 + (arr_azi(i)-1.5:0.03:arr_azi(i)+1.5-0.03);
    thetaRep = repmat(thetas,100,1);
    thetaRep1 = repmat(thetasActual,100,1);
    
    k_1 =  (4*pi/velLight*cosd(elevation)*fRep.*sind(thetaRep));
    k_2 =  (4*pi/velLight*cosd(elevation)*fRep.*cosd(thetaRep));

    k_1_2=  (4*pi/velLight*cosd(elevation)*fRep.*sind(thetaRep1));
    k_2_2 =  (4*pi/velLight*cosd(elevation)*fRep.*cosd(thetaRep1));
    mask1 = zeros(128,128);
    mask1(64-42:64+42,64-42:64+42) = 1;
    arr_img_comp(:,:,i)=img_comp;
    arr_img_comp_rot(:,:,i)=img_comp_rot;
    arr_img_fft(:,:,i) = fftshift(fft2(ifftshift(mask1.*img_comp)));
    centerIm = round(size(img_comp)/2);
    arr_img_fft_crop(:,:,i) =(1./taylorWindow.*arr_img_fft(centerIm(1)-50:centerIm(1)...
        + 49,centerIm(2)-49:centerIm(2) + 50 ,i));
    arr_img_fft_polar(:,:,i) = (interp2(XX,YY,(arr_img_fft_crop(:,:,i)),k_1,k_2,'nearest',0));
   
end

elevation=elevation*ones(size(arr_azi));
%Save data to disc
save(sprintf('../data/phase_histories/%s_PH',file_out_prefix),'arr_img_fft_polar','elevation','arr_azi','f','thetas');

