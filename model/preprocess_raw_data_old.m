% Load and match image with the given sample jpg
clear all;close all;clc;
% TODO: Specify path to the specific class' data from CD
fileNamePrefix = {'BMP2_SN_9563', 'BTR70_SN_C71', 'T72_SN_132', 'BTR_60','2S1', 'BRDM_2', 'D7', 'T62', 'ZIL131', 'ZSU_23_4'};
file_ext_arr={'000','004','015','003','000','001','005','016','025','026'};

for a=5:5%length(fileNamePrefix)
file_ext=file_ext_arr{a}
file_out_prefix = fileNamePrefix{a}%'T72_SN_132';  

pathLoad=sprintf('E:/BoxS/Box Sync/mstar/Data_from_CD/CD_2_publicTargets/TARGETS/TRAIN/17_DEG/%s/',file_out_prefix);

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
arr_img_fft_polar=c
%}

for i = 1:size(file_names,1)
    path1=[pathLoad file_names(i,:)];
    gg=MSTAR_LOAD_IMAGE(path1);
    depression = gg.MeasuredDepression;
%     yy = (4*pi/velLight*f*cosd(depression))*cosd(-1.5);
%     xx1 = (-1.5:0.03:1.5-0.03).';
%     xx = (4*pi/velLight*f(end)*cosd(depression)*sind(xx1));
%     [XX,YY] = meshgrid(xx,yy);
    
    img_comp=(flipud(gg.ImageData));
    size(img_comp)
    %figure;imagesc(abs(img_comp));
    
    azi=gg.TargetAz;
    %len=128;
    if i==1
    len=min(size(img_comp));
    end
    
    img_comp=img_comp(1:len,1:len);
    %figure;imagesc(abs(img_comp));
    img_comp_rot = imrotate(img_comp,azi,'bilinear','crop');

    arr_azi(i)=azi;
    
    thetasActual = 90 + (arr_azi(i)-1.5:0.03:arr_azi(i)+1.5-0.03);
    thetaRep = repmat(thetas,100,1);
    thetaRep1 = repmat(thetasActual,100,1);
    
    k_1 =  (4*pi/velLight*cosd(depression)*fRep.*sind(thetaRep));
    k_2 =  (4*pi/velLight*cosd(depression)*fRep.*cosd(thetaRep));
    yy = linspace(min(k_2(:)),max(k_2(:)),100);
    xx = linspace(min(k_1(:)),max(k_1(:)),100);
   
    [XX,YY] = meshgrid(xx,yy); 
    
    k_1_2=  (4*pi/velLight*cosd(depression)*fRep.*sind(thetaRep1));
    k_2_2 =  (4*pi/velLight*cosd(depression)*fRep.*cosd(thetaRep1));
    mask1 = zeros(len,len);
    mask1(64-31:64+32,64-31:64+32) = 1;
    arr_img_comp(:,:,i)=img_comp;
    arr_img_comp_rot(:,:,i)=img_comp_rot;
    
    %Masking then Coverting to polar domain
    arr_img_fft(:,:,i) = fftshift(fft2(ifftshift(mask1.*img_comp)));
    centerIm = round(size(img_comp)/2);
    %Undoing the Taylor window
    arr_img_fft_crop(:,:,i) =(1./taylorWindow.*arr_img_fft(centerIm(1)-49:centerIm(1)...
        + 50,centerIm(2)-49:centerIm(2) + 50 ,i));
    %figure; scatter(k_1(:),k_2(:));hold on;scatter(XX(:),YY(:),'r');
    % Cartesian to Polar
    arr_img_fft_polar(:,:,i) = (interp2(XX,YY,(arr_img_fft_crop(:,:,i)),k_1,k_2,'nearest',0));
    %arr_img_fft_polar(:,:,i) = griddata(XX,YY,real(arr_img_fft_crop(:,:,i)),k_1,k_2,'nearest')...
    %                        +1i*griddata(XX,YY,imag(arr_img_fft_crop(:,:,i)),k_1,k_2,'nearest');
    
end
end
%depression=depression*ones(size(arr_azi));
%Save data to disc
save(sprintf('../data/phase_histories/%s_PH',file_out_prefix),'arr_img_fft_polar','depression','arr_azi','f','thetas');