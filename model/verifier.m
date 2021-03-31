clc;
clear all;
close all;
fileNamePrefix = {'T72_SN_132'};
idxClass=1;
%% Load old data
PH_old=load(sprintf('../data/old/phase_histories/%s_masked',fileNamePrefix{idxClass}));
PH_new=load(sprintf('../data/phase_histories/%s_PH',fileNamePrefix{idxClass}));

C_old=load(sprintf('../data/old/recovered_coefficients/%s',fileNamePrefix{idxClass}));
C_new=load(sprintf('../data/recovered_coefficients/%s',fileNamePrefix{idxClass}));

%% Compare old vs new
isequal(PH_new.arr_img_fft_polar,PH_old.arr_img_fft_polar)
figure;imagesc(abs(PH_new.arr_img_fft_polar(:,:,1)));
figure;imagesc(abs(PH_old.arr_img_fft_polar(:,:,1)));

%% Compare trainIm vs. trainPH in the older SARopGen. Paths wrt mstar Box shared folder
PH_trainIm=load(sprintf('E:/BoxS/Box Sync/mstar/Tushar_tries/Augment_April_26_2020/train/%s_trainImPhaseHistories',fileNamePrefix{idxClass}));
PH_trainPH=load(sprintf('E:/BoxS/Box Sync/mstar/Tushar_tries/Augment_April_26_2020/train_phaseHistory/%s_masked',fileNamePrefix{idxClass}));

%% Load Final output
final_out=load(sprintf('../data/gen_aug_data/%s/sample_1',fileNamePrefix{idxClass}));
figure;imagesc(abs(permute(final_out.imgTrain(1,:,:),[3,2,1])));