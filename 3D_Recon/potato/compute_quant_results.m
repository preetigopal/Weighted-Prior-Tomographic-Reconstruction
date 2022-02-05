clear all;
clc;
close all;
outDirectory = 'results/TV_with_weighted_prior/18_views';

file = load(sprintf('%s/potato6_testVol.mat',outDirectory));
testVol = double(file.data);
file = load(sprintf('%s/FDK.mat',outDirectory));
fdkVol = double(file.FDK);
file = load(sprintf('%s/TV_lambdaTV_0.10.mat',outDirectory));
tvVol = file.TV;
file = load(sprintf('%s/weighted_prior_kk_1_lambda_prior_0.220000.mat',outDirectory));
weightedPriorVol = file.result_weighted_pca;


roi1_col = 55:55+15; % red RoI
roi1_row = 78:78+15;
roi2_col = 30:30+50;% green RoI
roi2_row = 60:60+35;
roi3_col = 20:20+80;
roi3_row = 40:40 +70;
% roi3_col = 130:145;
% roi3_row = 140:180;
roi1_frames = 30:37;
roi2_frames = 30:37;
roi3_frames = 30:37;


exponent_vector = [0.1 0.2 0.7];
 %exponent_vector = [1 1 1];

roi1_testVol = testVol(roi1_row,roi1_col,roi1_frames);
roi1_fbp = fdkVol(roi1_row,roi1_col,roi1_frames);
roi1_tv = tvVol(roi1_row,roi1_col,roi1_frames);
roi1_weighted_pca = weightedPriorVol(roi1_row,roi1_col,roi1_frames);


roi1_testVol = roi1_testVol - min(roi1_testVol(:));
roi1_testVol  = roi1_testVol./max(roi1_testVol(:));
roi1_fbp = roi1_fbp - min(roi1_fbp(:));
roi1_fbp = roi1_fbp./max(roi1_fbp(:));
roi1_tv = roi1_tv - min(roi1_tv(:));
roi1_tv = roi1_tv./max(roi1_tv(:));
roi1_weighted_pca = roi1_weighted_pca - min(roi1_weighted_pca(:));
roi1_weighted_pca = roi1_weighted_pca./max(roi1_weighted_pca(:));  

ssim1_fbp = ssim(roi1_fbp,roi1_testVol,'Exponent',exponent_vector);
ssim1_tv = ssim(roi1_tv,roi1_testVol,'Exponent',exponent_vector);
ssim1_weighted_pca = ssim(roi1_weighted_pca,roi1_testVol,'Exponent',exponent_vector); 

[rel_rmse1_fbp, rmse1_fbp] = rmseVol(roi1_fbp,roi1_testVol);
[rel_rmse1_tv, rmse1_tv] = rmseVol(roi1_tv,roi1_testVol);
[rel_rmse1_weighted_pca, rmse1_weighted_pca] = rmseVol(roi1_weighted_pca,roi1_testVol);

roi2_testVol = testVol(roi2_row,roi2_col,roi2_frames);
roi2_fbp = fdkVol(roi2_row,roi2_col,roi2_frames);
roi2_tv = tvVol(roi2_row,roi2_col,roi2_frames);
roi2_weighted_pca = weightedPriorVol(roi2_row,roi2_col,roi2_frames);

roi2_testVol = roi2_testVol - min(roi2_testVol(:));
roi2_testVol  = roi2_testVol./max(roi2_testVol(:));
roi2_fbp = roi2_fbp - min(roi2_fbp(:));
roi2_fbp = roi2_fbp./max(roi2_fbp(:));
roi2_tv = roi2_tv - min(roi2_tv(:));
roi2_tv = roi2_tv./max(roi2_tv(:));
roi2_weighted_pca = roi2_weighted_pca - min(roi2_weighted_pca(:));
roi2_weighted_pca = roi2_weighted_pca./max(roi2_weighted_pca(:));  

ssim2_fbp = ssim(roi2_fbp,roi2_testVol,'Exponent',exponent_vector);
ssim2_tv = ssim(roi2_tv,roi2_testVol,'Exponent',exponent_vector);
ssim2_weighted_pca = ssim(roi2_weighted_pca,roi2_testVol,'Exponent',exponent_vector); 

[rel_rmse2_fbp, rmse2_fbp] = rmseVol(roi2_fbp,roi2_testVol);
[rel_rmse2_tv, rmse2_tv] = rmseVol(roi2_tv,roi2_testVol);
[rel_rmse2_weighted_pca, rmse2_weighted_pca] = rmseVol(roi2_weighted_pca,roi2_testVol);

roi3_testVol = testVol(roi3_row,roi3_col,roi3_frames);
roi3_fbp = fdkVol(roi3_row,roi3_col,roi3_frames);
roi3_tv = tvVol(roi3_row,roi3_col,roi3_frames);
roi3_weighted_pca = weightedPriorVol(roi3_row,roi3_col,roi3_frames);

roi3_testVol = roi3_testVol - min(roi3_testVol(:));
roi3_testVol  = roi3_testVol./max(roi3_testVol(:));
roi3_fbp = roi3_fbp - min(roi3_fbp(:));
roi3_fbp = roi3_fbp./max(roi3_fbp(:));
roi3_tv = roi3_tv - min(roi3_tv(:));
roi3_tv = roi3_tv./max(roi3_tv(:));
roi3_weighted_pca = roi3_weighted_pca - min(roi3_weighted_pca(:));
roi3_weighted_pca = roi3_weighted_pca./max(roi3_weighted_pca(:));  

ssim3_fbp = ssim(roi3_fbp,roi3_testVol,'Exponent',exponent_vector);
ssim3_tv = ssim(roi3_tv,roi3_testVol,'Exponent',exponent_vector);
ssim3_weighted_pca = ssim(roi3_weighted_pca,roi3_testVol,'Exponent',exponent_vector); 

[rel_rmse3_fbp, rmse3_fbp] = rmseVol(roi3_fbp,roi3_testVol);
[rel_rmse3_tv, rmse3_tv] = rmseVol(roi3_tv,roi3_testVol);
[rel_rmse3_weighted_pca, rmse3_weighted_pca] = rmseVol(roi3_weighted_pca,roi3_testVol);
% 

% ssim_fbp = ssim(fdkVol(:,:,1:end-30),testVol(:,:,1:end-30),'Exponent',exponent_vector);
% ssim_tv = ssim(tvVol(:,:,1:end-30),testVol(:,:,1:end-30),'Exponent',exponent_vector);
% ssim_weighted_pca = ssim(weightedPriorVol(:,:,1:end-30),testVol(:,:,1:end-30),'Exponent',exponent_vector); 
% ssim_fbp = ssim(fdkVol(:,:,end),testVol(:,:,end));%,'Exponent',exponent_vector);
% ssim_tv = ssim(tvVol(:,:,end),testVol(:,:,end));%,'Exponent',exponent_vector);
% ssim_weighted_pca = ssim(weightedPriorVol(:,:,end),testVol(:,:,end));%,'Exponent',exponent_vector); 

testVol = testVol - min(testVol(:));
testVol  = testVol./max(testVol(:));
fdkVol = fdkVol - min(fdkVol(:));
fdkVol = fdkVol./max(fdkVol(:));
tvVol = tvVol - min(tvVol(:));
tvVol = tvVol./max(tvVol(:));
weightedPriorVol = weightedPriorVol - min(weightedPriorVol(:));
weightedPriorVol = weightedPriorVol./max(weightedPriorVol(:)); 

ssim_fbp = ssim(fdkVol,testVol,'Exponent',exponent_vector);
ssim_tv = ssim(tvVol,testVol,'Exponent',exponent_vector);
ssim_weighted_pca = ssim(weightedPriorVol,testVol,'Exponent',exponent_vector);

rmse_fbp = rmseVol(fdkVol,testVol);
rmse_tv = rmseVol(tvVol,testVol);
rmse_weighted_pca = rmseVol(weightedPriorVol,testVol);

sliceNo = 30;
testIm = testVol(:,:,sliceNo);
image = testIm;
image = image - min(image(:));
image = image./max(image(:));
fh = figure;imshow(image); hold on;
% rectangle('Position',[55, 78,15 ,15],'EdgeColor','r','LineWidth',2);hold on;
% rectangle('Position',[30, 60,50 ,35],'EdgeColor','g','LineWidth',2);%hold on;
rectangle('Position',[20, 40,80 ,70],'EdgeColor','r','LineWidth',2);
%name = sprintf('inputs/testIm.png');
name = sprintf('inputs/testIm_paper.png');
f=getframe; imwrite(f.cdata,name);





