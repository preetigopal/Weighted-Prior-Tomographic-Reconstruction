clear all;
clc;
close all;
% outDirectory = 'results/TV_with_weighted_prior/45_views_with_registered_templates';
outDirectory = 'results/TV_with_weighted_prior/30_views';
dataset ='okra3';

file = load(sprintf('%s/okra3_testVol.mat',outDirectory));
testVol = double(file.data);
% testVol = testVol(:,110:end-100,:);
file = load(sprintf('%s/FDK.mat',outDirectory));
fdkVol = double(file.FDK);
% fdkVol = fdkVol(:,110:end-100,:);
file = load(sprintf('%s/TV_lambdaTV_0.45.mat',outDirectory));
tvVol = file.TV;
% tvVol = tvVol(:,110:end-100,:);
file = load(sprintf('%s/weighted_prior_kk_1_lambda_prior_0.500000.mat',outDirectory));
weightedPriorVol = file.result_weighted_pca;
% weightedPriorVol = weightedPriorVol(:,110:end-100,:);





% testVol = testVol - min(testVol(:));
% testVol  = testVol./max(testVol(:));
% fdkVol = fdkVol - min(fdkVol(:));
% fdkVol = fdkVol./max(fdkVol(:));
% tvVol = tvVol - min(tvVol(:));
% tvVol = tvVol./max(tvVol(:));
% weightedPriorVol = weightedPriorVol - min(weightedPriorVol(:));
% weightedPriorVol = weightedPriorVol./max(weightedPriorVol(:));  


% leaveSlices = 0;
% ssim_fdk = ssim(fdkVol(:,:,end-leaveSlices),testVol(:,:,end-leaveSlices))
% ssim_tv = ssim(tvVol(:,:,end-leaveSlices),testVol(:,:,end-leaveSlices))
% ssim_wightedPrior = ssim(weightedPriorVol(:,:,end-leaveSlices),testVol(:,:,end-leaveSlices))
% 
% relRMSE_fdk = rmseVol(fdkVol,testVol)
% relRMSE_tv = rmseVol(tvVol,testVol)
% relRMSE_weightedPrior = rmseVol(weightedPriorVol,testVol)
% 
% ssim_fdk_1 = ssim(fdkVol(:,:,end-leaveSlices),testVol(:,:,end-leaveSlices),'Exponent',[0.1 0.2 0.7])
% ssim_tv_1 = ssim(tvVol(:,:,end-leaveSlices),testVol(:,:,end-leaveSlices),'Exponent',[0.1 0.2 0.7])
% ssim_wightedPrior_1 = ssim(weightedPriorVol(:,:,end-leaveSlices),testVol(:,:,end-leaveSlices),'Exponent',[0.1 0.2 0.7])


roi1_col = 145:160;
roi1_row = 171:197;
roi2_col = 160:200;%160:185
roi2_row = 131:131+86;
% roi3_col = 130:145;
% roi3_row = 140:180;
roi1_frames = 30:37;
roi2_frames = 32:38;

ssim_exponent = [0.1 0.2 0.7];
%ssim_exponent = [1 1 1];


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

ssim1_fbp = ssim(roi1_fbp,roi1_testVol,'Exponent',ssim_exponent);
ssim1_tv = ssim(roi1_tv,roi1_testVol,'Exponent',ssim_exponent);
ssim1_weighted_pca = ssim(roi1_weighted_pca,roi1_testVol,'Exponent',ssim_exponent); 
[relRMSE1_fbp, rmse1_fbp] = rmseVol(roi1_fbp,roi1_testVol)
[relRMSE1_tv, rmse1_tv] = rmseVol(roi1_tv,roi1_testVol)
[relRMSE1_weighted_pca, rmse1_weighted_pca] = rmseVol(roi1_weighted_pca,roi1_testVol)


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

ssim2_fbp = ssim(roi2_fbp,roi2_testVol,'Exponent',ssim_exponent);
ssim2_tv = ssim(roi2_tv,roi2_testVol,'Exponent',ssim_exponent);
ssim2_weighted_pca = ssim(roi2_weighted_pca,roi2_testVol,'Exponent',ssim_exponent); 
[relRMSE2_fbp, rmse2_fbp] = rmseVol(roi2_fbp,roi2_testVol)
[relRMSE2_tv, rmse2_tv] = rmseVol(roi2_tv,roi2_testVol)
[relRMSE2_weighted_pca, rmse2_weighted_pca] = rmseVol(roi2_weighted_pca,roi2_testVol)


testVol = testVol - min(testVol(:));
testVol  = testVol./max(testVol(:));
fdkVol = fdkVol - min(fdkVol(:));
fdkVol = fdkVol./max(fdkVol(:));
tvVol = tvVol - min(tvVol(:));
tvVol = tvVol./max(tvVol(:));
weightedPriorVol = weightedPriorVol - min(weightedPriorVol(:));
weightedPriorVol = weightedPriorVol./max(weightedPriorVol(:));  

% relRMSE_fbp = rmseVol(fdkVol,testVol);
% relRMSE_tv = rmseVol(tvVol,testVol);
% relRMSE_weighted_pca = rmseVol(weightedPriorVol,testVol);

ssim_fbp = ssim(fdkVol,testVol,'Exponent',ssim_exponent);
ssim_tv = ssim(tvVol,testVol,'Exponent',ssim_exponent);
ssim_weighted_pca = ssim(weightedPriorVol,testVol,'Exponent',ssim_exponent); 


sliceNo = 31;
testIm = testVol(:,:,sliceNo);


image = testIm;
image = image - min(image(:));
image = image./max(image(:));
image = image(:,110:end-100);
fh = figure;imshow(image); hold on;
rectangle('Position',[145-110, 171,16 ,20],'EdgeColor','r','LineWidth',2);%hold on;
%rectangle('Position',[160-110, 131,40 ,86],'EdgeColor','g','LineWidth',2);
name = sprintf('inputs/testImCropped.png');

f=getframe; imwrite(f.cdata,name);
