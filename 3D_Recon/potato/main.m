% potato_few_views

dataset = 'potato6';
outDirectory = ('results/TV_with_weighted_prior/18_views/');
mkdir(outDirectory);
projFolderName =  sprintf('inputs/%s/projf32_LIN_AP_nc',dataset);
dirInfo = dir(projFolderName);
sliceNumber = 30;
lambda_prior_list = [0.1 0.2 0.22 0.25 0.3 0.4 1];

%% -------------Read Projection data---------------------------

for j = 3:size(dirInfo,1)
    fileName = dirInfo(j).name;
    filePath = sprintf('%s/%s',projFolderName,fileName);
    ncid = netcdf.open(filePath);
    projections = double(netcdf.getVar(ncid,1));
    if j==3
        proj = projections;
    else
        proj = cat(3,proj,projections);                  
    end
    netcdf.close(ncid);
end

  
startz = 1; endz = size(proj,3); % corresponds to each view
starty = 121; endy = size(proj,2)-121+1;
startx = 125; endx = size(proj,1)-125+1; 

proj = proj(startx:endx, starty: endy, startz:endz); 
%selectedProj = 1:20:size(proj,3);
selectedProj = 1:50:size(proj,3);
proj = proj(:,:,selectedProj);
proj1 = zeros(150,150,size(proj,3));
for i = 1:size(proj,3)
    image = imresize(proj(:,:,i),[150 150]);
    proj1(:,:,i) = image;
end
proj = proj1;
angStepSize = 20;
numAngles = length(selectedProj);
ParamSetting_potato_pca(angStepSize,outDirectory)
name = sprintf('%s/parameters_potato_pca.mat',outDirectory);
load(name);

dim = [param.nx param.ny param.nz]; 

sizeProj = size(proj);
save('ProjectionSize.mat','sizeProj');


%% --------generate Ground truth for final comparison------------------

name = sprintf('inputs/%s_cropped_900views_fdk',dataset);
testVolFile = load(name);
testVol = testVolFile.FDK;
name = 'testVol';
store_video(testVol,name,outDirectory,dataset);

name = sprintf('%s/testVol.png',outDirectory);
temp = testVol(:,:,sliceNumber);
temp = temp - min(temp(:));
temp = temp/max(temp(:));
imwrite(temp,name);

testVol1 = testVol - min(testVol(:));
testVol1  = testVol1./max(testVol1(:));



%% -----Filtered Back-projection ---------------------------------------

disp('starting FDK');
tic;
[FDK] = potato_FDK(proj,param,dim);
time_FDK = toc;

name = 'FDK';
store_video(FDK,name,outDirectory,dataset);
name = sprintf('%s/FDK.mat',outDirectory);
save(name,'testVol','FDK');
name = sprintf('%s/FDK.png',outDirectory);
temp = FDK(:,:,sliceNumber);
temp = temp - min(temp(:));
temp = temp/max(temp(:));
imwrite(temp,name);

FDK1 = FDK - min(FDK(:));
FDK1 = FDK1./max(FDK1(:));


%% -----------TV (No object-prior) ---------------------------------

m = size(proj(:),1);
A_both_functions =@(z,mode) systemFunction(z,mode,dim,param,m);
y = reshape(proj,[sizeProj(1)*sizeProj(2)*sizeProj(3) 1]);



lambdaTV_list = [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30];%0.10;%0.05:0.05:1.1;
for i = 1:length(lambdaTV_list)
    lambdaTV = lambdaTV_list(i);
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

%% ------------- TV + Weighted object-prior--------------------------
          
disp('Starting weighted prior reconstructions');
baseMethod = 'TV';
lambda_for_base_method = 0;%lambdaTV_list(1);
lambda_for_pilot_recon = {1,0.10};%{1,lambdaTV_list(1)};
kk = 1;
templateNos = [1,3,4];
genEigenSpacePotato_hole(templateNos,outDirectory)

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








