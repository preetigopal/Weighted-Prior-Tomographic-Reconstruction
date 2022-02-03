function driver(dataset)
%%
% This code is part of the following work which has been submitted to IEEE Transactions on Computational Imaging for peer review.

% Title: "Eliminating object prior-bias from sparse-projection tomographic reconstructions"
% Authors: Preeti Gopal, Sharat Chandran, Imants Svalbe and Ajit Rajwade 

%  This is the main file. Run this with the arguement 'tmh_7' for
%  reconstructing 2D Liver.


if strcmp(dataset,'potato')
    
    angleSet =  [12]    ;
    name = sprintf('../../../3D/potato/inputs/potato%d_cropped_900views_fdk.mat',6) ;
    fileData = load(name);
    volume = fileData.FDK;
    testIm= double(volume(:,:,30));    
    lambdaTV_list = 0.05;%:0.05:1.1;
    lambda_prior_list = [0.1 1.1 2.1 3.1 4.1 5.1]; % The actual lambda_prior value is squared of this value.
    kk_values = 30;
    lambda_TV = 0.05;  
    
elseif strcmp(dataset,'okra')
    
    sliceNum = 30;
    angleSet = [48];    
    name = sprintf('../../data/templates/okra/okra6_okra%d_reg_450views_fdk.mat',3) ;
    fileData = load(name);
    volume = fileData.FDK;
    testIm= double(volume(15:end-14,15:end-14,sliceNum));    
    lambdaTV_list = 0.01;%:0.005:0.1;
    lambda_prior_list = 0.7;%[0.01, 0.03, 0.05, 0.5, 0.7, 0.9, 1.1, 5, 10]; % The actual lambda_prior value is squared of this value.
    kk_values = [5,10,20,40,50,60,70,80,90,100];%,370,410,450];%30;
    lambda_TV = 0.01;
    
elseif strcmp(dataset,'sprouts')
    
    sliceNum = 38;
    angleSet = 20;%[20, 30];    
    name = sprintf('../../3D/Sprouts/inputs/sprouts%d_givenVolume.mat',11) ;
    fileData = load(name);
    volume = fileData.volume;
    testIm= double(volume(10:end-4,2:end-2,38));
%     testIm= double(volume(:,:,sliceNum)); 
    lambdaTV_list = 0.1;%:0.05:0.30
    lambda_prior_list = 2.1;%[0.1 1.1 2.1 3.1 4.1]; % The actual lambda_prior value is squared of this value.
    kk_values = 30;
    lambda_TV = 0.1;
    
elseif strcmp(dataset,'tmh_7')
    
    angleSet = 130;% [30];
    testNo = 16  ;
    dataFileName = 'RFA2 - Anon/Rfa/RFA 1.5 B30s -'
    sliceNum = 27;
 
    projFolderName =  sprintf('%s %d/',dataFileName,testNo);
    dirInfo = dir(projFolderName);
    for j = 1:size(dirInfo,1)-2
        fileName = dirInfo(j+2).name;
        filePath = sprintf('%s/%s',projFolderName,fileName);
        image = dicomread(filePath);
        if j==1
            testVol = zeros([size(image,1) size(image,2) (size(dirInfo,1)-2)]);
        end
        testVol(:,:,j) = image;
    end
    testIm= testVol(:,:,sliceNum);
    testIm = testIm - min(testIm(:));
    testIm = testIm./max(testIm(:));
    
    lambdaTV_list = 0.036;
    lambda_prior_list = 1.1;
    kk_values = 30;
    lambda_TV = 0.036;  
    
elseif strcmp(dataset,'tmh_8')
    
    angleSet = [60]%[120];
    testNo = 18  ;
    dataFileName = 'RFA2 - Anon/Rfa/RFA 1.5 B30s -'
    sliceNum = 27;
 
    projFolderName =  sprintf('../tmh/%s %d/',dataFileName,testNo);
    dirInfo = dir(projFolderName);
    for j = 1:size(dirInfo,1)-2
        fileName = dirInfo(j+2).name;
        filePath = sprintf('%s/%s',projFolderName,fileName);
        image = dicomread(filePath);
        if j==1
            testVol = zeros([size(image,1) size(image,2) (size(dirInfo,1)-2)]);
        end
        testVol(:,:,j) = image;
    end
    testIm= testVol(:,:,sliceNum);
    testIm = testIm - min(testIm(:));
    testIm = testIm./max(testIm(:));
    
    lambdaTV_list = 0.038;
    lambda_prior_list = [0.05 0.1 0.15 1.1 2.1 3.1]; 
    kk_values = 30;
    lambda_TV = 0.038; 
    
end
%%
for iter=1:length(angleSet)
    numAngles = angleSet(iter);
    outDirectory = sprintf('results/%s/%d_views/slice_%d/',dataset,numAngles,sliceNum);
    mkdir(outDirectory);
    numAngles = angleSet(iter);
    [y,idx1,dim] = generateMeasurements(testIm,numAngles,outDirectory);
    
    %Pilot Reconstruction Methods
    % 0: FBP
    % 1: TV
    disp('Starting pilot reconstructions');
    pilotReconMethods = [0,1]; 
    lambda_for_pilot_recon = {1,lambdaTV_list};
    performPilotReconstruction(y,idx1,dim,pilotReconMethods,lambda_for_pilot_recon,outDirectory);
    
    disp('Starting weighted prior reconstructions');
    baseMethod = 'TV';
    lambda_for_base_method = lambda_TV;
    lambda_for_pilot_recon = {1,lambdaTV_list(1)};
    for kk_iter = 1:length(kk_values)
        kk = kk_values(kk_iter)
        performWeightedPriorReconstruction(dataset,y,idx1,dim,pilotReconMethods,lambda_for_pilot_recon, baseMethod,...
            lambda_for_base_method,lambda_prior_list,kk,sliceNum,outDirectory)
    end
end

end

