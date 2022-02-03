clear all;
close all;

output_dir = 'results/TV_with_weighted_prior/45_views_with_registered_templates';
name = sprintf('%s/okra3_testVol.mat',output_dir);
f1 = load(name);

name = sprintf('%s/FDK.mat',output_dir);
f3 = load(name);

name = sprintf('%s/TV_lambdaTV_0.50.mat',output_dir);
f4 = load(name);

name = sprintf('%s/weighted_prior_kk_1_lambda_prior_0.400000.mat',output_dir);
f6 = load(name);



%-------------store the images-------------

fdkIm = f1.data(:,:,30);
ImCropped = fdkIm(:,110:end-100);
ImCropped = ImCropped - min(ImCropped(:));
ImCropped = ImCropped./max(ImCropped(:));
imwrite(ImCropped,'test_cropped.png');


fdkIm = f3.FDK(:,:,30);
ImCropped = fdkIm(:,110:end-100);
ImCropped = ImCropped - min(ImCropped(:));
ImCropped = ImCropped./max(ImCropped(:));
imwrite(ImCropped,'fdk_cropped.png');


tvIm = f4.TV(:,:,30);
ImCropped = tvIm(:,110:end-100);
ImCropped = ImCropped - min(ImCropped(:));
ImCropped = ImCropped./max(ImCropped(:));
imwrite(ImCropped,'TV_cropped.png');


prior_weightedIm = f6.result_weighted_pca(:,:,30);
ImCropped = prior_weightedIm(:,110:end-100);
ImCropped = ImCropped - min(ImCropped(:));
ImCropped = ImCropped./max(ImCropped(:));
imwrite(ImCropped,'prior_weighted_cropped.png');

