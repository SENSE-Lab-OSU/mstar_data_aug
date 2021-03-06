clear all;close all;clc;
%stuff from main code
fileNamePrefix = {'2S1','BMP2_SN_9563','BMP2_SN_9566','BMP2_SN_C21',...
                'BRDM_2','BTR_60','BTR70_SN_C71','D7','T62',...
                'T72_SN_132','T72_SN_812','T72_SN_S7','ZIL131','ZSU_23_4'};
%n_samples=[233,233,232,256,299,298,299,299,299,299]
array_dtheta=[-6:0.25:-0.25 0.25:0.25:6];
shiftsy =[0,0.15,0, 0.15];
shiftsx =[0,0, 0.15,0.15];


idxClass=1; %specify class here
path2coeff=('../data/recovered_coefficients/');
m=load(sprintf('%s%s',path2coeff,fileNamePrefix{idxClass}));
numTrainingSamples = size(m.y_recovered,3);
clearvars m

%Merge code
infer_sample_flag=0;
sample_start=1;
sample_end=numTrainingSamples;
class_folder=sprintf('../data/gen_aug_data/%s',fileNamePrefix{idxClass});
class_file=sprintf('../data/gen_aug_data/%s_aug_images.mat',fileNamePrefix{idxClass});
s=load(sprintf('%s/sample_%d',class_folder,sample_start));
factr=size(s.elev,1);

%%
%Check if partial file exists and load it
if isfile(class_file)
    load(class_file);
    %find apparent sample_start
    sample_start_hat=find(aziTrain==0, 1 ); %first idx that is zero
    sample_start_hat=1+((sample_start_hat-1)/factr);
    if rem(sample_start_hat,1)~=0
        error('Sample start from file MUST be an integer but is %d.',...
            sample_start_hat);
    end
    
    if infer_sample_flag==1
        sample_start=sample_start_hat;
    end
    
else
    %Initialize empty array
    imgTrain =zeros((length(array_dtheta)+1)*length(shiftsy)*numTrainingSamples,64,64);
    aziTrain=zeros((length(array_dtheta)+1)*length(shiftsy)*numTrainingSamples,1);
    elev=zeros((length(array_dtheta)+1)*length(shiftsy)*numTrainingSamples,1);
    if sample_start~=1
        error('Sample starts do not match and are: specified= %d, from_file= %d.',...
            sample_start,1);
    end
end
%%
%fill the complete matrices
fprintf('Starting filling class %d from sample %d\n',idxClass,sample_start);
for idxSample=sample_start:sample_end
    fill_start=1+(idxSample-1)*factr;
    fill_end=fill_start+(factr-1);
    fprintf('Filling sample %d in array indices %d to %d\n',idxSample,fill_start,fill_end);
    %load and fill
    s=load(sprintf('%s/sample_%d',class_folder,idxSample));
    imgTrain(fill_start:fill_end,:,:)=s.imgTrain;
    aziTrain(fill_start:fill_end,:)=s.aziTrain;
    elev(fill_start:fill_end,:)=s.elev;
end

%% Save
disp('Saving Data to disk...')
save(class_file,'imgTrain','aziTrain','elev','-v7.3');