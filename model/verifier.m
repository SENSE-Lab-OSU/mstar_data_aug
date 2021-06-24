clc;
clear all;
close all;
fileNamePrefix = {'2S1','BMP2_SN_9563','BMP2_SN_9566','BMP2_SN_C21',...
                'BRDM_2','BTR_60','BTR70_SN_C71','D7','T62',...
                'T72_SN_132','T72_SN_812','T72_SN_S7','ZIL131','ZSU_23_4'};
idxClass=1;
%figure;imagesc(abs(squeeze(imgTrain(20,:,:))));
%% Load old data
%PH_old=load(sprintf('../data/old/phase_histories/%s_masked',fileNamePrefix{idxClass}));
%PH_old=load(sprintf('E:/BoxS/Box Sync/mstar/Tushar_tries/Augment_April_26_2020/train/%s_trainImPhaseHistories',fileNamePrefix{idxClass}));
PH_old=load(sprintf('E:/BoxS/Box Sync/mstar/Tushar_tries/Augment_April_26_2020/train_phaseHistory/%s_masked',fileNamePrefix{idxClass}));
PH_new=load(sprintf('../data/phase_histories/%s_PH',fileNamePrefix{idxClass}));

C_old=load(sprintf('../data/old/recovered_coefficients/%s',fileNamePrefix{idxClass}));
C_new=load(sprintf('../data/recovered_coefficients/%s',fileNamePrefix{idxClass}));

%% Compare old vs new
close all;
isequal(PH_new.arr_azi,PH_old.arr_azi)
isequal(PH_new.depression,PH_old.depression)
isequal(PH_new.arr_img_fft_polar,PH_old.arr_img_fft_polar)
%a1=C_new.x_recovered(:,2:end)==0;prod(a1(:))

img_idx=1;
figure;imagesc(abs(PH_new.arr_img_fft_polar(:,:,img_idx)));title('PH_{new}.arr-img-fft-polar')
figure;imagesc(abs(PH_old.arr_img_fft_polar(:,:,img_idx)));title('PH_{old}.arr-img-fft-polar')

diff_arr=abs(PH_old.arr_img_fft_polar(:,:,img_idx)-PH_new.arr_img_fft_polar(:,:,img_idx))/abs(PH_new.arr_img_fft_polar(:,:,img_idx));
disp(min(abs(diff_arr(:))));disp(mean(abs(diff_arr(:))));disp(max(abs(diff_arr(:))));

figure;imagesc(diff_arr);title('perc diff_{arr}');
figure;imagesc(abs(squeeze(PH_new.arr_img_comp(img_idx,:,:))));title('True base image');

test_imag_ifft=ifftshift(ifft2(PH_old.arr_img_fft_polar(:,:,img_idx)));
test_imag_ifft_new=ifftshift(ifft2(PH_new.arr_img_fft_polar(:,:,img_idx)));

test_imag_ifft=test_imag_ifft./norm(test_imag_ifft,'fro');
test_imag_ifft_new=test_imag_ifft_new./norm(test_imag_ifft_new,'fro');

figure;imagesc(abs(squeeze((test_imag_ifft_new(:,:)))));title('normed ifft(PH_{new})');
figure;imagesc(abs(squeeze((test_imag_ifft(:,:)))));title('normed ifft(PH_{old})');

diff_ifft=test_imag_ifft_new-test_imag_ifft;
disp(min(abs(diff_ifft(:))));disp(mean(abs(diff_ifft(:))));disp(max(abs(diff_ifft(:))));

figure;imagesc(abs(squeeze((diff_ifft))));title('normed ifft(PH_{new})-ifft(PH_{old})');

%test_imag_ifft=ifftshift(ifft2(PH_old.arr_img_fft_polar(:,:,1)));
%figure;imagesc(abs(squeeze(flipud(test_imag_ifft(:,:)))));

%% Compare recovered coefficients
%Q about RC size
ax_min=min([abs(C_old.x_recovered(:,img_idx));abs(C_new.x_recovered(:,img_idx))]);
ax_max=max([abs(C_old.x_recovered(:,img_idx));abs(C_new.x_recovered(:,img_idx))]);

figure;
plot(abs(C_old.x_recovered(:,img_idx)),...
    abs(C_new.x_recovered(:,img_idx)),'ro',...
    [ax_min,ax_max],[ax_min,ax_max],'k','linewidth',2);
grid on;title('Comparison of old C (x-axis) vs. new (y-axis) C');
xlabel('old C');ylabel('new C');
%plot(C_new.x_recovered(:,img_idx),'r-+');legend('Old C','New C');
%% Compare trainIm vs. trainPH in the older SARopGen. Paths wrt mstar Box shared folder
PH_trainIm=load(sprintf('E:/BoxS/Box Sync/mstar/Tushar_tries/Augment_April_26_2020/train/%s_trainImPhaseHistories',fileNamePrefix{idxClass}));
PH_trainPH=load(sprintf('E:/BoxS/Box Sync/mstar/Tushar_tries/Augment_April_26_2020/train_phaseHistory/%s_masked',fileNamePrefix{idxClass}));

%% Load Final output
final_out=load(sprintf('../data/gen_aug_data/%s/sample_1',fileNamePrefix{idxClass}));
%mask1 = zeros(numPixelsCrop,numPixelsCrop);
centerIm=64;
aoi=final_out.aziTrain(1);

figure;
im_preproc=abs(squeeze(PH_new.arr_img_comp(img_idx,centerIm-32:centerIm+31,centerIm-32:centerIm+31)));
im_preproc=im_preproc./norm(im_preproc,'fro');
imagesc(20*log10(im_preproc));colormap('jet');
title(sprintf('Cropped True base image at %s deg.',...
    num2str(round(aoi,2))));grid on;
caxis_lims=[max(20*log10(im_preproc(:)))-30 max(20*log10(im_preproc(:)))];
caxis(caxis_lims);

%figure;imagesc(squeeze(abs(final_out.imgTrain(1,:,:))));title('algorithm output at ');
for idx_aug=15:3:49
figure;
azi_shift=aoi-final_out.aziTrain(idx_aug);
im_preproc=squeeze(abs(final_out.imgTrain(idx_aug,:,:)));
im_preproc=im_preproc./norm(im_preproc,'fro');
imagesc(20*log10(im_preproc));colormap('jet');
title(sprintf('Augmented Image at %s deg. azi_shifted %s deg.',...
    num2str(round(aoi,2)),num2str(round(azi_shift,2))));grid on;
%caxis([max(20*log10(im_preproc(:)))-30 max(20*log10(im_preproc(:)))]);
caxis(caxis_lims);
end