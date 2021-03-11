Contains the main code to produce synthetic complex SAR images for data augmentation.

You need MATLAB (preferably 2017+) and python 3 currently to proceed. Feel free to contribute and improve compatibility.


1. Run "python data/get_data.py". This should download required data into two folders, "phase_histories" and "recovered_coefficients" needed for now to proceed. (We are working on uploading complete end to end code, i.e. the code to generate the files in these folders using publicly available MSTAR data from the CD.)

Further Image Generation is done class by class to enable execution on multiple PCs.

2. Open "model/generate_aug_images.m" and set some parameters as per need. Set variable "idxClass" to the class number you want to generate for, where class nos. are in reference to the variable "fileNamePrefix". Use "sample_start" variable to start generating from that sample and check comments to use "sample_end" similarily. Run this script. This should start generating images with subpixel shifts (refer paper) image by image for that particular class and saving it in "data/gen_aug_data/<class-name>".

Each sample's mat file will have [1 (original sample)+ 48 (extrapolated samples)] * 4(subpixel shifts)]= 196 complex-valued images in -1.5 to 1.5 deg neighborhood of that sample's azimuthal angle. For reference, azimuthal angles are provided in the variable "aziTrain" correpsonding to the complex images in "imgTrain" where the first image is the original.

3. Once done, use the "model/merge_files.m" to combine all the individual samples into a single array for easier data-loading later. Set variables "idxClass", "sample_start" and "sample_end" just like in step 2. Run this script. This should start merging all samples for that particular class into a final single array and saving it as "data/gen_aug_data/<class-name>_aug_images.mat". Once done, you may delete the corresponding directory "data/gen_aug_data/<class-name>".