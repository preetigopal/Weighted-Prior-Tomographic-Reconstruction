

function [eigenVecs,meanTemplate] = genEigenSpace2Dtmh(templateNos,sliceNum,dataFileName,numSlices,sizeOfSlice,outDirectory)


% This code is part of the following work which has been submitted to IEEE Transactions on Computational Imaging for peer review.

% Title: "Eliminating object prior-bias from sparse-projection tomographic reconstructions"
% Authors: Preeti Gopal, Sharat Chandran, Imants Svalbe and Ajit Rajwade 

% This function computes high quality eigen-space from the prior-images
% of the liver.

numTemplates = length(templateNos);
volumes  = zeros(sizeOfSlice(1),sizeOfSlice(2),numSlices,numTemplates);

% Read the templates in 3D format-------------------

for i = 1:numTemplates  
   %name = sprintf('../../data/templates/potato/potato4_potato%d_reg_900views_fdk.mat',templateNos(i)) ;
   
   %projFolderName =  sprintf('%s %d/',dataFileName,templateNos(i));
   projFolderName =  sprintf('%s %d/',dataFileName,templateNos(i));
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
      volumes(:,:,:,i) = testVol;
end

% Get the templates and compute their mean----------------------------------------------------------
templateIm = cell(1,numTemplates);
for i = 1:numTemplates
    
    in = volumes(:,:,sliceNum,i); 
    in = in - min(in(:));
    in  = in./max(in(:));
    temp = in; % double(in);

    

    %imshow(temp,[]); title(i); pause(0.5);    
    if i==1
        templates = zeros((size(temp,1)*size(temp,2)), numTemplates);
        templates(:,1) = temp(:);
        sumI = zeros(size(templates(:,i)));
        
        tt = temp;
        minimum = min(tt(:));
        tt = tt - minimum;
        maximum = max(tt(:));
    end
    templates(:,i) = temp(:);
    sumI = sumI + templates(:,i);
    
    fname = sprintf('%s/template%d.mat',outDirectory,i);
    save(fname,'temp');
    
    temp = temp - minimum;
    temp = temp./maximum;
    templateIm{i} = temp;
    
    name = sprintf('%s/template%d.png',outDirectory,i);
    temp = templateIm{i};
    imwrite(temp,name); 
end

meanTemplate = sumI./numTemplates;
%figure;imshow(reshape(meanTemplate,size(temp)),[]);title('Mean image');


% Compute the Covariance Matrix--------------------------------------------------------------------

templates = templates -repmat(meanTemplate,1,numTemplates);
L = templates'*templates;
[W,D] = eig(L);
V = templates*W;
V = normc_fcn(V);
[m n] = size(V);

% picking top k eigen values and their corresponding vectors-----------------------------------------------------
% This forms the eigen space of the covariance matrix of the templates-----------------                  

numDim = numTemplates-1;
eigenVals = zeros(1,numDim);
eigenVecs = zeros(m,numDim);
%figure;
for j = 1:numDim    
    eigenVals(j) = D(n-j+1,n-j+1);
    eigenVecs(:,j) = V(:,n-j+1);
 %   imshow(reshape(eigenVecs(:,j),size(temp)),[]);title(j);pause(1);
end

name = sprintf('%s/EigenSpace_2D_tmh.mat',outDirectory);
save(name,'eigenVals','eigenVecs','meanTemplate','minimum','maximum');

%-------------------------------------------------------------------------
%-------------Testing quality of Eigen space constructed---------------------------
%------------------------------------------------------------------------
% Compute the weights ('alpha') for each of the templates (you're taking the
% projection (dot product) of each template onto each eigenvector
% Note: Here- cols of eigenVecs are eigenvectors.

alpha = cell(1,numTemplates);

for i = 1:numTemplates
    alpha{i} = eigenVecs'*(templates(:,i));  
end

% Reconstruct each template back from its projection coefficients on the
% eigen vectors

for i = 1:numTemplates
    coeff = alpha{i};
    recon = zeros(size(meanTemplate));
    for j = 1:numDim
        recon = recon + (coeff(j)*eigenVecs(:,j));
    end
    recon = recon + meanTemplate;
   % imshow([reshape(templates(:,i)+meanTemplate,size(temp)) reshape(recon,size(temp))],[]);title(i);pause(.5)
end

