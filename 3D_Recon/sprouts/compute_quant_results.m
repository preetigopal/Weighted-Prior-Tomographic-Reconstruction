clear all;
clc;
close all;
outDirectory = 'results/TV_with_weighted_prior/30_views';


file = load(sprintf('%s/sprouts11_testVol.mat',outDirectory));
testVol = double(file.data);
file = load(sprintf('%s/FDK.mat',outDirectory));
fdkVol = double(file.FDK);
file = load(sprintf('%s/TV_lambdaTV_0.15.mat',outDirectory));
tvVol = file.TV;
file = load(sprintf('%s/weighted_prior_kk_1_lambda_prior_0.500000.mat',outDirectory));
weightedPriorVol = file.result_weighted_pca;

roi1_col = 90:115;
roi1_row = 75:95;
roi1_frames = 33:39;

exponent_factor = [0.1 0.2 0.7]
%  exponent_factor = [1 1 1]

roi1_testVol = testVol(roi1_row,roi1_col,roi1_frames);
roi1_fbp = fdkVol(roi1_row,roi1_col,roi1_frames);
roi1_tv = tvVol(roi1_row,roi1_col,roi1_frames);
roi1_weighted_pca = weightedPriorVol(roi1_row,roi1_col,roi1_frames);

roi1_testVol  = roi1_testVol  - min(roi1_testVol (:));
roi1_testVol  = roi1_testVol ./max(roi1_testVol (:));
roi1_fbp = roi1_fbp - min(roi1_fbp(:));
roi1_fbp = roi1_fbp./max(roi1_fbp(:));
roi1_tv = roi1_tv - min(roi1_tv(:));
roi1_tv = roi1_tv./max(roi1_tv(:));
roi1_weighted_pca = roi1_weighted_pca - min(roi1_weighted_pca(:));
roi1_weighted_pca = roi1_weighted_pca./max(roi1_weighted_pca(:));  

ssim1_fbp = ssim(roi1_fbp,roi1_testVol,'Exponent',exponent_factor);
ssim1_tv = ssim(roi1_tv,roi1_testVol,'Exponent',exponent_factor);
ssim1_weighted_pca = ssim(roi1_weighted_pca,roi1_testVol,'Exponent',exponent_factor); 
[relRMSE1_fbp, rmse1_fbp] = rmseVol(roi1_fbp,roi1_testVol)
[relRMSE1_tv, rmse1_tv] = rmseVol(roi1_tv,roi1_testVol)
[relRMSE1_weighted_pca, rmse1_weighted_pca] = rmseVol(roi1_weighted_pca,roi1_testVol)


testVol = testVol - min(testVol(:));
testVol  = testVol./max(testVol(:));
fdkVol = fdkVol - min(fdkVol(:));
fdkVol = fdkVol./max(fdkVol(:));
tvVol = tvVol - min(tvVol(:));
tvVol = tvVol./max(tvVol(:));
weightedPriorVol = weightedPriorVol - min(weightedPriorVol(:));
weightedPriorVol = weightedPriorVol./max(weightedPriorVol(:)); 


ssim_fbp = ssim(fdkVol,testVol,'Exponent',exponent_factor);
ssim_tv = ssim(tvVol,testVol,'Exponent',exponent_factor);
ssim_weighted_pca = ssim(weightedPriorVol,testVol,'Exponent',exponent_factor);    


%------------saving one of the result slices

sliceNo = 38;
testIm = testVol(:,:,sliceNo);

image = testIm;
fh = figure;imshow(image,[]); hold on;
rectangle('Position',[90, 75,26 ,21],'EdgeColor','r','LineWidth',2);hold on;
name = sprintf('inputs/testIm_paper.png');
f=getframe; imwrite(f.cdata,name);








