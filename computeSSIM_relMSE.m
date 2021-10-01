clc;
clear all; close all;
dataset = 'sprouts'

if strcmp(dataset,'tmh_7')
    
    outDirectory = 'results/tmh_7/30_views/';
    fileName = sprintf('%s/testIm.mat',outDirectory);
    data = load(fileName);
    testIm = data.testIm;
    
    roi_row = 135:135+140;
    roi_col = 275:275+155;
    testIm_roi = testIm(roi_col,roi_row);
    testIm_roi =  testIm_roi - min(testIm_roi(:));
    testIm_roi =  testIm_roi./max(testIm_roi(:));
    
    
    testIm_zoomed = testIm(275:275+155,135:135+140);    
    testIm_zoomed  = testIm_zoomed  - min(testIm_zoomed (:));
    testIm_zoomed = testIm_zoomed ./max(testIm_zoomed (:));
    name = sprintf('%s/test_zoomed.png',outDirectory);
    imwrite(testIm_zoomed,name);
%     figure;imshow(testIm_roi,[]);
%     fh = figure;imshow(testIm,[]);
%     impixelinfo;
%     hold on;
%     rectangle('Position',[135, 275, 140, 155],'EdgeColor','g','LineWidth',3);
%     name = sprintf('%s/testIm_roi.png',outDirectory);
%     f=getframe; imwrite(f.cdata,name);

    fileName = sprintf('%s/fbp.mat',outDirectory);
    data = load(fileName);
    fbp = data.result_fbp;
    fbp_roi = fbp(roi_col,roi_row);
    fbp_roi =  fbp_roi - min(fbp_roi(:));
    fbp_roi =  fbp_roi./max(fbp_roi(:));
    
    fbp_zoomed = fbp(275:275+155,135:135+140);       
    fbp_zoomed  = fbp_zoomed  - min(fbp_zoomed (:));
    fbp_zoomed = fbp_zoomed ./max(fbp_zoomed (:));
    name = sprintf('%s/fbp_zoomed.png',outDirectory);
    imwrite(fbp_zoomed,name);
    
    relMSE_FBP =  mean2((fbp_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
    ssimVal_FBP = ssim(fbp_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    
%     lambdaTV_list1 = [0.010,0.020,0.030,0.040,0.050,0.060,0.070,0.080];
%     relMSE1 = zeros(1,length(lambdaTV_list1));
%     ssimVal1 = zeros(1,length(lambdaTV_list1));
    
    lambdaTV_list2 = 0.036;%[0.032,0.034,0.036,0.038];
    relMSE2 = zeros(1,length(lambdaTV_list2));
    ssimVal2 = zeros(1,length(lambdaTV_list2));
    
%     for i=1:length(lambdaTV_list1)
%         lambdaTV = lambdaTV_list1(i)
%         fileName = sprintf('%s/tv_lambdaTV_%.3f.mat',outDirectory,lambdaTV);
%         data = load(fileName);
%         result = data.result_tv;
%         result_roi = result(roi_col,roi_row);
%         relMSE1(1,i) =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
%         ssimVal1(1,i) = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
%     end
    
    for i=1:length(lambdaTV_list2)
        lambdaTV = lambdaTV_list2(i)
        fileName = sprintf('%s/tv_lambdaTV_%.3f.mat',outDirectory,lambdaTV);
        data = load(fileName);
        result = data.result_tv;
        result_roi = result(roi_col,roi_row);
        result_roi = result_roi - min(result_roi(:));
        result_roi = result_roi./max(result_roi(:));
        
        tv_zoomed = result(275:275+155,135:135+140);       
        tv_zoomed  = tv_zoomed  - min(tv_zoomed (:));
        tv_zoomed = tv_zoomed ./max(tv_zoomed (:));
        name = sprintf('%s/tv_zoomed.png',outDirectory);
        imwrite(tv_zoomed,name);
        
        relMSE2(1,i) =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
        ssimVal2(1,i) = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    end   
%     
    lambdaPrior_list = 1.1;%[0.1,1.1,2.1,3.1,4.1];
    relMSE_prior = zeros(1,length(lambdaPrior_list));
    ssimVal_prior = zeros(1,length(lambdaPrior_list));
    

    for i=1:length(lambdaPrior_list)
        lambda_prior = lambdaPrior_list(i)
        fileName = sprintf('%s/weighted_prior_kk_30_lambda_prior_%f.mat',outDirectory,lambda_prior);
        data = load(fileName);
        result = data.result_weighted_pca;
        result_roi = result(roi_col,roi_row);
        result_roi = result_roi - min(result_roi(:));
        result_roi = result_roi./max(result_roi(:));   
        
        weighted_pca_zoomed = result(275:275+155,135:135+140);       
        weighted_pca_zoomed  = weighted_pca_zoomed  - min(weighted_pca_zoomed (:));
        weighted_pca_zoomed = weighted_pca_zoomed ./max(weighted_pca_zoomed (:));
        name = sprintf('%s/weighted_pca_zoomed.png',outDirectory);
        imwrite(weighted_pca_zoomed,name);
        
        relMSE_prior(1,i) =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
        ssimVal_prior(1,i) = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    end
    
    
end
%---------------------------------------------------------------------------------------------------------------------
if strcmp(dataset,'tmh_8')
    
    outDirectory = 'results/tmh_8/60_views/';
    fileName = sprintf('%s/testIm.mat',outDirectory);
    data = load(fileName);
    testIm = data.testIm;
    
    roi_row = 135:135+140;
    roi_col = 275:275+155;
    testIm_roi = testIm(roi_col,roi_row);
    testIm_roi =  testIm_roi - min(testIm_roi(:));
    testIm_roi =  testIm_roi./max(testIm_roi(:));
    figure;imshow(testIm_roi,[]);
    fh = figure;imshow(testIm,[]);
%     impixelinfo;
%     hold on;
%     rectangle('Position',[135, 275, 140, 155],'EdgeColor','g','LineWidth',3);
%     name = sprintf('%s/testIm_roi.png',outDirectory);
%     f=getframe; imwrite(f.cdata,name);

    fileName = sprintf('%s/fbp.mat',outDirectory);
    data = load(fileName);
    fbp = data.result_fbp;
    fbp_roi = fbp(roi_col,roi_row);
    fbp_roi =  fbp_roi - min(fbp_roi(:));
    fbp_roi =  fbp_roi./max(fbp_roi(:));
    relMSE_FBP =  mean2((fbp_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
    ssimVal_FBP = ssim(fbp_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    
    lambdaTV_list = [0.038];%[0.010,0.020,0.030,0.032,0.034,0.036,0.038,0.040,0.050,0.060,0.070,0.080];
    relMSE = zeros(1,length(lambdaTV_list));
    ssimVal = zeros(1,length(lambdaTV_list));
    
    
    for i=1:length(lambdaTV_list)
        lambdaTV = lambdaTV_list(i)
        fileName = sprintf('%s/tv_lambdaTV_%.3f.mat',outDirectory,lambdaTV);
        data = load(fileName);
        result = data.result_tv;
        result_roi = result(roi_col,roi_row);
        result_roi = result_roi - min(result_roi(:));
        result_roi = result_roi./max(result_roi(:));
        relMSE(1,i) =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
        ssimVal(1,i) = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    end
%     
%         lambdaPrior_list = 1.1;%[0.1,1.1,2.1,3.1,4.1];
%     relMSE_prior = zeros(1,length(lambdaPrior_list));
%     ssimVal_prior = zeros(1,length(lambdaPrior_list));
    
    lambdaPrior_list = [0.05,0.1,0.15,1.1,2.1];
    relMSE_prior = zeros(1,length(lambdaPrior_list));
    ssimVal_prior = zeros(1,length(lambdaPrior_list));

    for i=1:length(lambdaPrior_list)
        lambda_prior = lambdaPrior_list(i)
        fileName = sprintf('%s/weighted_prior_kk_30_lambda_prior_%f.mat',outDirectory,lambda_prior);
        data = load(fileName);
        result = data.result_weighted_pca;
        result_roi = result(roi_col,roi_row);
        result_roi = result_roi - min(result_roi(:));
        result_roi = result_roi./max(result_roi(:));
        relMSE_prior(1,i) =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
        ssimVal_prior(1,i) = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    end
    
 
    
    
end

if strcmp(dataset,'potato')
    
    roi_col = 20:20+80;
    roi_row = 40:40 +70;
    
    outDirectory = 'results/potato/6_views/';
    fileName = sprintf('%s/testIm.mat',outDirectory);
    data = load(fileName);
    testIm = data.testIm;
    
    testIm_roi = testIm(roi_row,roi_col);
    testIm_roi =  testIm_roi - min(testIm_roi(:));
    testIm_roi =  testIm_roi./max(testIm_roi(:));
    figure;imshow(testIm_roi,[]);
    
    testIm_zoomed = testIm(40:40 +70,20:20+80);    
    testIm_zoomed  = testIm_zoomed  - min(testIm_zoomed (:));
    testIm_zoomed = testIm_zoomed ./max(testIm_zoomed (:));
    name = sprintf('%s/test_zoomed.png',outDirectory);
    imwrite(testIm_zoomed,name);
    
    fileName = sprintf('%s/fbp.mat',outDirectory);
    data = load(fileName);
    fbp = data.result_fbp;
    fbp_roi = fbp(roi_row,roi_col);
    fbp_roi =  fbp_roi - min(fbp_roi(:));
    fbp_roi =  fbp_roi./max(fbp_roi(:));
    
    fbp_zoomed = fbp(40:40 +70,20:20+80);   
    fbp_zoomed  = fbp_zoomed  - min(fbp_zoomed (:));
    fbp_zoomed = fbp_zoomed ./max(fbp_zoomed (:));
    name = sprintf('%s/fbp_zoomed.png',outDirectory);
    imwrite(fbp_zoomed,name);
    
    relMSE_FBP =  mean2((fbp_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
    ssimVal_FBP = ssim(fbp_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);

    lambdaTV = 0.05
    fileName = sprintf('%s/tv_lambdaTV_%.2f.mat',outDirectory,lambdaTV);
    data = load(fileName);
    result = data.result_tv;
    result_roi = result(roi_row,roi_col);
    result_roi = result_roi - min(result_roi(:));
    result_roi = result_roi./max(result_roi(:));
    
    tv_zoomed = result(40:40 +70,20:20+80);   
    tv_zoomed  = tv_zoomed  - min(tv_zoomed (:));
    tv_zoomed = tv_zoomed ./max(tv_zoomed (:));
    name = sprintf('%s/tv_zoomed.png',outDirectory);
    imwrite(tv_zoomed,name);
    
    relMSE_TV =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
    ssimVal_TV = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    
    lambda_prior = 1.1
    fileName = sprintf('%s/weighted_prior_kk_30_lambda_prior_%f.mat',outDirectory,lambda_prior);
    data = load(fileName);
    result = data.result_weighted_pca;
    result_roi = result(roi_row,roi_col);
    result_roi = result_roi - min(result_roi(:));
    result_roi = result_roi./max(result_roi(:));
    
    weighted_pca_zoomed = result(40:40 +70,20:20+80);   
    weighted_pca_zoomed  = weighted_pca_zoomed  - min(weighted_pca_zoomed (:));
    weighted_pca_zoomed = weighted_pca_zoomed ./max(weighted_pca_zoomed (:));
    name = sprintf('%s/weighted_pca_zoomed.png',outDirectory);
    imwrite(weighted_pca_zoomed,name);
        
    relMSE_prior =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
    ssimVal_prior = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    
end


if strcmp(dataset,'okra')
    
    outDirectory = 'results/okra/48_views/slice_30/';
    fileName = sprintf('%s/testIm.mat',outDirectory);
    data = load(fileName);
    testIm = data.testIm;
%     
%     roi_row = 120:120+50;
%     roi_col = 120:120+80;

    roi_col = 142-9:160-9;
    roi_row = 160:190;
    
    testIm_roi = testIm(roi_row,roi_col);
    figure;imshow(testIm_roi,[]);
    testIm_roi =  testIm_roi - min(testIm_roi(:));
    testIm_roi =  testIm_roi./max(testIm_roi(:));
    
    testIm_zoomed = testIm(150:190,113:171);    
    testIm_zoomed  = testIm_zoomed  - min(testIm_zoomed (:));
    testIm_zoomed = testIm_zoomed ./max(testIm_zoomed (:));
    name = sprintf('%s/test_zoomed.png',outDirectory);
    imwrite(testIm_zoomed,name);
%     figure;imshow(testIm_roi,[]);
%     fh = figure;imshow(testIm,[]);
%     impixelinfo;
%     hold on;
% %     rectangle('Position',[120, 120, 50, 80],'EdgeColor','g','LineWidth',3);
%     name = sprintf('%s/testIm_roi.png',outDirectory);
%     f=getframe; imwrite(f.cdata,name);
    

    fileName = sprintf('%s/fbp.mat',outDirectory);
    data = load(fileName);
    fbp = data.result_fbp;
    fbp_roi = fbp(roi_row,roi_col);
    fbp_roi =  fbp_roi - min(fbp_roi(:));
    fbp_roi =  fbp_roi./max(fbp_roi(:));
    
    fbp_zoomed = fbp(150:190,113:171);    
    fbp_zoomed  = fbp_zoomed  - min(fbp_zoomed (:));
    fbp_zoomed = fbp_zoomed ./max(fbp_zoomed (:));
    name = sprintf('%s/fbp_zoomed.png',outDirectory);
    imwrite(fbp_zoomed,name);
    
    relMSE_FBP =  mean2((fbp_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
    ssimVal_FBP = ssim(fbp_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);

    lambdaTV = 0.01
    fileName = sprintf('%s/tv_lambdaTV_%.3f.mat',outDirectory,lambdaTV);
    data = load(fileName);
    result = data.result_tv;
    result_roi = result(roi_row,roi_col);
    result_roi = result_roi - min(result_roi(:));
    result_roi = result_roi./max(result_roi(:));
    
    tv_zoomed = result(150:190,113:171);    
    tv_zoomed  = tv_zoomed  - min(tv_zoomed (:));
    tv_zoomed = tv_zoomed ./max(tv_zoomed (:));
    name = sprintf('%s/tv_zoomed.png',outDirectory);
    imwrite(tv_zoomed,name);
    
    relMSE_TV =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
    ssimVal_TV = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    
    lambdaPrior_list = 0.7;%[0.5,0.7,0.9,1.1];
    relMSE_prior = zeros(1,length(lambdaPrior_list));
    ssimVal_prior = zeros(1,length(lambdaPrior_list));

    for i=1:length(lambdaPrior_list)
        lambda_prior = lambdaPrior_list(i)
        fileName = sprintf('%s/weighted_prior_kk_30_lambda_prior_%f.mat',outDirectory,lambda_prior);
        data = load(fileName);
        result = data.result_weighted_pca;
        result_roi = result(roi_row,roi_col);
        result_roi = result_roi - min(result_roi(:));
        result_roi = result_roi./max(result_roi(:));
        
        weighted_pca_zoomed = result(150:190,113:171);    
        weighted_pca_zoomed  = weighted_pca_zoomed  - min(weighted_pca_zoomed (:));
        weighted_pca_zoomed = weighted_pca_zoomed ./max(weighted_pca_zoomed (:));
        name = sprintf('%s/weighted_pca_zoomed.png',outDirectory);
        imwrite(weighted_pca_zoomed,name);
        
        relMSE_prior(1,i) =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
        ssimVal_prior(1,i) = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    end

%     lambda_prior = 0.7;
%     kk_list = [5,10,20,30,40,50,60,70,80,90,100,130];
%     relMSE_kk = zeros(1,length(kk_list));
%     ssimVal_kk = zeros(1,length(kk_list));
%     
%     for i=1:length(kk_list)
%         kk = kk_list(i)
%         fileName = sprintf('%s/weighted_prior_kk_%d_lambda_prior_%f.mat',outDirectory,kk,lambda_prior);
%         data = load(fileName);
%         result = data.result_weighted_pca;
%         result_roi = result(roi_row,roi_col);
%         result_roi = result_roi - min(result_roi(:));
%         result_roi = result_roi./max(result_roi(:));
%         
%         weighted_pca_zoomed = result(150:190,113:171);    
%         weighted_pca_zoomed  = weighted_pca_zoomed  - min(weighted_pca_zoomed (:));
%         weighted_pca_zoomed = weighted_pca_zoomed ./max(weighted_pca_zoomed (:));
%         name = sprintf('%s/weighted_pca_zoomed.png',outDirectory);
%         imwrite(weighted_pca_zoomed,name);
%         
%         relMSE_kk(1,i) =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
%         ssimVal_kk(1,i) = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
%     end
%     
end


if strcmp(dataset,'sprouts')
    
    %outDirectory = 'results/sprouts/20_views/slice_38/';
    outDirectory = 'results/sprouts/100_views/';
    fileName = sprintf('%s/testIm.mat',outDirectory);
    data = load(fileName);
    testIm = data.testIm;
    
    roi_col = 180:220;
    roi_row = 155:190;
    
    testIm_roi = testIm(roi_row,roi_col);
    testIm_roi = testIm_roi - min(testIm_roi(:));
    testIm_roi = testIm_roi./max(testIm_roi(:));
    figure;imshow(testIm_roi,[]);
    
    testIm_zoomed = testIm(105:195,160:240);    
    testIm_zoomed  = testIm_zoomed  - min(testIm_zoomed (:));
    testIm_zoomed = testIm_zoomed ./max(testIm_zoomed (:));
    name = sprintf('%s/test_zoomed.png',outDirectory);
    imwrite(testIm_zoomed,name);
    
    fileName = sprintf('%s/fbp.mat',outDirectory);
    data = load(fileName);
    fbp = data.result_fbp;
    fbp_roi = fbp(roi_row,roi_col);
    fbp_roi = fbp_roi - min(fbp_roi(:));
    fbp_roi = fbp_roi./max(fbp_roi(:));
    
    fbp_zoomed = fbp(105:195,160:240);    
    fbp_zoomed  = fbp_zoomed  - min(fbp_zoomed (:));
    fbp_zoomed = fbp_zoomed ./max(fbp_zoomed (:));
    name = sprintf('%s/fbp_zoomed.png',outDirectory);
    imwrite(fbp_zoomed,name);
    
    relMSE_FBP =  mean2((fbp_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
    ssimVal_FBP = ssim(fbp_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);

    lambdaTV = 0.1
    fileName = sprintf('%s/tv_lambdaTV_%.3f.mat',outDirectory,lambdaTV);
    data = load(fileName);
    result = data.result_tv;
    result_roi = result(roi_row,roi_col);
    result_roi = result_roi - min(result_roi(:));
    result_roi = result_roi./max(result_roi(:));
    
    tv_zoomed = result(105:195,160:240);    
    tv_zoomed  = tv_zoomed  - min(tv_zoomed (:));
    tv_zoomed = tv_zoomed ./max(tv_zoomed (:));
    name = sprintf('%s/tv_zoomed.png',outDirectory);
    imwrite(tv_zoomed,name);
    
    relMSE_TV =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
    ssimVal_TV = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    ssimVal_TV_standard = ssim(result_roi,testIm_roi);
    
    lambdaPrior_list = [2.1]%[0.1,1.1,2.1,3.1,4.1];
    relMSE_prior = zeros(1,length(lambdaPrior_list));
    ssimVal_prior = zeros(1,length(lambdaPrior_list));

    for i=1:length(lambdaPrior_list)
        lambda_prior = lambdaPrior_list(i)
        fileName = sprintf('%s/weighted_prior_kk_30_lambda_prior_%f.mat',outDirectory,lambda_prior);
        data = load(fileName);
        result = data.result_weighted_pca;
        result_roi = result(roi_row,roi_col);
        result_roi = result_roi - min(result_roi(:));
        result_roi = result_roi./max(result_roi(:));
        
        weighted_pca_zoomed = result(105:195,160:240);    
        weighted_pca_zoomed  = weighted_pca_zoomed  - min(weighted_pca_zoomed (:));
        weighted_pca_zoomed = weighted_pca_zoomed ./max(weighted_pca_zoomed (:));
        name = sprintf('%s/weighted_pca_zoomed.png',outDirectory);
        imwrite(weighted_pca_zoomed,name);
        
        relMSE_prior(1,i) =  mean2((result_roi - testIm_roi).^2)/mean2(testIm_roi.^2);
        ssimVal_prior(1,i) = ssim(result_roi,testIm_roi,'Exponents',[0.1,0.2,0.7]);
    end
    


    image = testIm;
    fh = figure;imshow(image,[]); hold on;
    % rectangle('Position',[90, 75,26 ,21],'EdgeColor','r','LineWidth',2);hold on;
    % rectangle('Position',[160-110, 131,40 ,86],'EdgeColor','g','LineWidth',2);
    rectangle('Position',[180, 155,41 ,36],'EdgeColor','r','LineWidth',2);
    name = sprintf('%s/testIm_roi.png',outDirectory);
    f=getframe; imwrite(f.cdata,name);   

end

if strcmp(dataset,'potato_all_methods')
    numAngles = 6;
    outDirectory = sprintf('results/comparing_all_methods/%d_angles',numAngles);
    fileName = sprintf('%s/testIm.mat',outDirectory);
    data = load(fileName);
    testIm = data.testIm;
    
    roi_col = 1:150;%20:20+80;
    roi_row = 1:150;%40:40 +70;
    
    roi_testIm = testIm(roi_row,roi_col);
    roi_testIm = roi_testIm - min(roi_testIm(:));
    roi_testIm = roi_testIm./max(roi_testIm(:));
        
    name = sprintf('%s/fbp.mat',outDirectory);
    data = load(name);
    fbp = data.result_fbp;
    
    roi_fbp = fbp(roi_row,roi_col);
    roi_fbp = roi_fbp - min(roi_fbp(:));
    roi_fbp = roi_fbp./max(roi_fbp(:));
    
    name = sprintf('%s/tv.mat',outDirectory);
    data = load(name);
    tv = data.result_tv;
    
    roi_tv = tv(roi_row,roi_col);
    roi_tv = roi_tv - min(roi_tv(:));
    roi_tv = roi_tv./max(roi_tv(:));
    
    name = sprintf('%s/cs_wavelet.mat',outDirectory);
    data = load(name);
    cs_wavelet = data.result_cs_wavelet;
    
    roi_cs_wavelet = cs_wavelet(roi_row,roi_col);
    roi_cs_wavelet = roi_cs_wavelet - min(roi_cs_wavelet(:));
    roi_cs_wavelet = roi_cs_wavelet./max(roi_cs_wavelet(:));
    
    name = sprintf('%s/cs_dct.mat',outDirectory);
    data = load(name);
    cs_dct = data.result_cs_dct;
    
    roi_cs_dct = cs_dct(roi_row,roi_col);
    roi_cs_dct = roi_cs_dct - min(roi_cs_dct(:));
    roi_cs_dct = roi_cs_dct./max(roi_cs_dct(:));
    
    name = sprintf('%s/art.mat',outDirectory);
    data = load(name);
    art = data.result_art;
    
    roi_art = art(roi_row,roi_col);
    roi_art = roi_art - min(roi_art(:));
    roi_art = roi_art./max(roi_art(:));
    
    name = sprintf('%s/sart.mat',outDirectory);
    data = load(name);
    sart = data.result_sart;
    
    roi_sart = sart(roi_row,roi_col);
    roi_sart = roi_sart - min(roi_sart(:));
    roi_sart = roi_sart./max(roi_sart(:));
    
    name = sprintf('%s/sirt.mat',outDirectory);
    data = load(name);
    sirt = data.result_sirt;
    
    roi_sirt = sirt(roi_row,roi_col);
    roi_sirt = roi_sirt - min(roi_sirt(:));
    roi_sirt = roi_sirt./max(roi_sirt(:));
    
    kk = 30;
    name = sprintf('%s/weighted_pca_all_methods%d.mat',outDirectory,kk);
    data = load(name);
    all_methods = data.result_weighted_pca;
    
    roi_all_methods = all_methods(roi_row,roi_col);
    roi_all_methods = roi_all_methods - min(roi_all_methods(:));
    roi_all_methods = roi_all_methods./max(roi_all_methods(:));
    
    figure;imshow([testIm, tv, all_methods],[]);
    figure;imshow([testIm all_methods],[]);
    
    ssimWts = [0.1,0.2,0.7];
    
    ssim_fbp = ssim(roi_fbp,roi_testIm,'Exponent',ssimWts)
    ssim_cs_dct = ssim(roi_cs_dct,roi_testIm,'Exponent',ssimWts)
    ssim_cs_wavelet = ssim(roi_cs_wavelet,roi_testIm,'Exponent',ssimWts)
    ssim_art = ssim(roi_art,roi_testIm,'Exponent',ssimWts)
    ssim_sart = ssim(roi_sart,roi_testIm,'Exponent',ssimWts)
    ssim_sirt = ssim(roi_sirt,roi_testIm,'Exponent',ssimWts)
    ssim_tv = ssim(roi_tv,roi_testIm,'Exponent',ssimWts)
    % ssim_plain_pca = ssim(roi_plain_pca,roi_testIm,'Exponent',ssimWts);
    ssim_all_methods = ssim(roi_all_methods,roi_testIm,'Exponent',ssimWts) 
    
end











