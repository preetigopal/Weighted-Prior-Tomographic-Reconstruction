close all;
clear all;
clc;


outDirectory = 'results/TV_with_weighted_prior/30_views';
dataset = 'okra3';

file = load('inputs/okra3_cropped_450views_fdk.mat');
testVol = file.FDK;
testVol = testVol(:,:,1:60);
testVol = testVol - min(testVol(:));
testVol = testVol./max(testVol(:));
name = sprintf('%s/FDK.mat',outDirectory);
file = load(name);
fdkVol = file.FDK;
fdkVol = fdkVol(:,:,1:60);
fdkVol = fdkVol - min(fdkVol(:));
fdkVol = fdkVol./max(fdkVol(:));
name = sprintf('%s/TV_lambdaTV_0.45.mat',outDirectory);
file = load(name);
tvVol = file.TV;
tvVol = tvVol(:,:,1:60);
tvVol = tvVol - min(tvVol(:));
tvVol = tvVol./max(tvVol(:));
name = sprintf('%s/weighted_prior_kk_1_lambda_prior_0.500000.mat',outDirectory);
file = load(name);
weightedPriorVol = file.result_weighted_pca;
weightedPriorVol = weightedPriorVol(:,:,1:60);
weightedPriorVol = weightedPriorVol - min(weightedPriorVol(:));
weightedPriorVol = weightedPriorVol./max(weightedPriorVol(:));



name = 'all_results';
all_results = [testVol fdkVol tvVol weightedPriorVol];
store_video(all_results,name,outDirectory,dataset);
    
