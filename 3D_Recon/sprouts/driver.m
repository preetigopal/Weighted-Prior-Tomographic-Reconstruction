clc;clear all;

dataset = 'sprouts11'
outDirectory = ('results/TV_with_weighted_prior/30_views/');
mkdir(outDirectory);
sliceNumber = 38;
lambda_prior_list = 0.5;


angStepSize = 12;
ParamSetting_sprouts(angStepSize)
load('parameters_sprouts.mat');
dim = [param.nx param.ny param.nz]; 



%% Generate synthetic projections------------------------

name = 'inputs/sprouts_11_syntheticVolume_1800views.mat';
file = load(name);
volume = double(file.data);
proj = double(CTprojection(volume,param));
sizeProj = size(proj);
save('ProjectionSize.mat','sizeProj');

%% --------generate Ground truth------------------

testVol = volume;
name = 'testVol';
store_video(testVol,name,outDirectory,dataset);

name = sprintf('%s/testVol.png',outDirectory);
temp = testVol(:,:,sliceNumber);
temp = temp - min(temp(:));
temp = temp./max(temp(:));
imwrite(temp,name);

testVol1 = testVol - min(testVol(:));
testVol1  = testVol1./max(testVol1(:));



%% -----Filtered Back-projection-----------------------------------------------------------------

disp('starting FDK');
tic;
[FDK] = sprouts_FDK(proj,param,dim);
time_FDK = toc;

name = 'FDK';
store_video(FDK,name,outDirectory,dataset);
name = sprintf('%s/FDK.mat',outDirectory);
save(name,'testVol','FDK');
name = sprintf('%s/FDK.png',outDirectory);
temp = FDK(:,:,sliceNumber);
temp = temp - min(temp(:));
temp = temp./max(temp(:));
imwrite(temp,name);

FDK1 = FDK - min(FDK(:));
FDK1 = FDK1./max(FDK1(:));


%% --------TV Reconstruction (No object-prior)
m = size(proj(:),1);
A_both_functions =@(z,mode) systemFunction(z,mode,dim,param,m);
y = reshape(proj,[sizeProj(1)*sizeProj(2)*sizeProj(3) 1]);



lambdaTV_list = 0.17;
for i = 1:length(lambdaTV_list)
    lambdaTV = lambdaTV_list(i)
    pilotReconMethods = [0,1]; 
    
    disp('starting TV');    
    fname = sprintf('%s/TV_lambdaTV_%.2f.mat',outDirectory,lambdaTV);            
    if(isfile(fname))
        datafile = load(fname);
        TV = datafile.TV;
        TV1 = TV - min(TV(:));
        TV1 = TV1./max(TV1(:));
        continue;
    else
        TV = plain_TV(y,dim,lambdaTV,A_both_functions);
        save(fname,'TV');
        name = sprintf('TV_lambdaTV_%.2f',lambdaTV); 
        store_video(TV,name,outDirectory,dataset);
        name = sprintf('%s/TV_lambdaTV_%.2f.png',outDirectory,lambdaTV);
        temp = TV(:,:,sliceNumber);
        temp = temp - min(temp(:));
        temp = temp/max(temp(:));
        imwrite(temp,name);
        
        TV1 = TV - min(TV(:));
        TV1 = TV1./max(TV1(:));
        all_volumes = [testVol1 FDK1 TV1];
        summary_video_name = sprintf('all_volumes_no_prior_lambdaTV_%.2f',lambdaTV); 
        store_video(all_volumes,summary_video_name,outDirectory,dataset);
    end
end


%% -----------------TV + weighted object-prior reconstruction----------------------
          
disp('Starting weighted prior reconstructions');
baseMethod = 'TV';
lambda_for_base_method = lambdaTV_list(1)
lambda_for_pilot_recon = {1,lambdaTV_list(1)};
kk = 1;
templateNos = [6,7,8,9,10];
genEigenSpaceSprouts_hole(templateNos,outDirectory)

for iter = 1:length(lambda_prior_list)
    lambda_prior = lambda_prior_list(iter)
    [result_weighted_pca,weightsIm] = performWeightedPriorReconstruction(dataset,y,param,dim,pilotReconMethods,lambda_for_pilot_recon,...
        baseMethod,lambda_for_base_method,lambda_prior,kk,sliceNumber,outDirectory);
    fname = sprintf('%s/weighted_prior_kk_%d_lambda_prior_%f.mat',outDirectory,kk,lambda_prior);
    save(fname,'result_weighted_pca');
    name = sprintf('%s/weighted_prior_kk_%d_lambda_prior_%f.png',outDirectory,kk,lambda_prior);
    temp = result_weighted_pca(:,:,sliceNumber);
    temp = temp - min(temp(:));
    temp = temp./max(temp(:));
    imwrite(temp,name);    

    fname = sprintf('%s/weightsIm_kk_%d_lambda_prior_%f.mat',outDirectory,kk,lambda_prior);
    save(fname,'weightsIm');
    name = sprintf('%s/weightsIm_kk_%d_lambda_prior_%f.png',outDirectory,kk,lambda_prior);
    temp = weightsIm(:,:,sliceNumber);
    temp = temp - min(temp(:));
    temp = temp./max(temp(:));
    imwrite(temp,name);  
    
    result_weighted_pca1 = result_weighted_pca - min(result_weighted_pca(:));
    result_weighted_pca1 = result_weighted_pca1./max(result_weighted_pca1(:));
    all_volumes = [testVol1 FDK1 TV1 result_weighted_pca1];
    summary_video_name = sprintf('all_volumes_weighted_prior_lambda_prior_%.2f',lambda_prior); 
    store_video(all_volumes,summary_video_name,outDirectory,dataset);
end








