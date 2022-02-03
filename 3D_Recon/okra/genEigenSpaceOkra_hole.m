function genEigenSpaceOkra_hole(templateNos,outDirectory)

startTemplateNum = templateNos(1);
endTemplateNum = templateNos(end);
numTemplates = length(startTemplateNum:endTemplateNum);
numSlices = 123; 
sizeOfSlice = [338 338];

volumes  = zeros(sizeOfSlice(1),sizeOfSlice(2),numSlices,numTemplates);

% Read the templates in 3D format-------------------

for i = 1:numTemplates  
   name = sprintf('inputs/okra%d_okra6_reg_cropped_450views_fdk.mat',i+startTemplateNum-1) 
   % --The templates have been registered with repect to Okra6 (6th volume of Okra)

   fileData = load(name);
   volumes(:,:,:,i) = fileData.FDK;
   image = volumes(:,:,30,i);
   image = image - min(image(:));
   image = image./max(image(:));
   name = sprintf('inputs/template_%d.png',i);
   imwrite(image,name);
   
   % store cropped image
   
   image = volumes(:,110:end-100,30,i);
   image = image - min(image(:));
   image = image./max(image(:));
   name = sprintf('inputs/templateCropped_%d.png',i);
   imwrite(image,name);
end


% store templates and test together------------------

name = sprintf('inputs/okra%d_cropped_450views_fdk.mat',3) ;
fileDataTest = load(name);
testVol = fileDataTest.FDK;
image = testVol(:,110:end-100,30);
image = image - min(image(:));
image = image./max(image(:));
name = 'inputs/test.png';
imwrite(image,name);



all_templates = [volumes(:,:,1:end-65,1) volumes(:,:,1:end-65,2) volumes(:,:,1:end-65,3) volumes(:,:,1:end-65,4) testVol(:,:,1:end-65)];
all_templates = all_templates - min(all_templates(:));
all_templates = all_templates./max(all_templates(:));
name = 'templates_and_test';
dataset ='_';
store_video(all_templates,name,'inputs',dataset);
tempTestImg = [volumes(:,:,10,1) volumes(:,:,10,2) volumes(:,:,10,3) volumes(:,:,10,4) testVol(:,:,10)];
tempTestImg = tempTestImg - min(tempTestImg(:));
tempTestImg = tempTestImg./max(tempTestImg(:));
imwrite(tempTestImg,'inputs/template_test_image.png');



% resize the templates--------------------------------------------------------
templates = zeros(sizeOfSlice(1)*sizeOfSlice(2)*numSlices,numTemplates);

for i = 1:numTemplates  
   temp = volumes(:,:,:,i);
   temp = temp(:);   
   
   templates(:,i) = temp;
end

% Get the templates and compute their mean----------------------------------------

for i = 1:numTemplates
    if i==1
        sumI = zeros(size(templates(:,i)));        
        tt = temp;
        minimum = min(tt(:));
        tt = tt - minimum;
        maximum = max(tt(:));
    end
    sumI = sumI + templates(:,i);
end
meanTemplate = sumI./numTemplates;


% Compute the Covariance Matrix--------------------------------------------------


templates = templates -repmat(meanTemplate,1,numTemplates);
L = templates'*templates;
[W,D] = eig(L);
V = templates*W;
V = normc(V);
[m n] = size(V);

% picking top k eigen values and their corresponding vectors-----------------------
% This forms the eigen space of the covariance matrix of the templates-----------------                  

numDim = numTemplates;%-1;
eigenVals = zeros(1,numDim);
eigenVecs = zeros(m,numDim);

for j = 1:numDim    
    eigenVals(j) = D(n-j+1,n-j+1);
    eigenVecs(:,j) = V(:,n-j+1);
    eigenImage = reshape(eigenVecs(:,j),[sizeOfSlice(1) sizeOfSlice(2) numSlices]);
%     for k = 1:numSlices
%         figure;imshow(eigenImage(:,:,k),[]);
%     end
end
name = sprintf('%s/EigenSpace_okra_hole.mat',outDirectory);
save(name,'eigenVals','eigenVecs','meanTemplate','minimum','maximum','-v7.3');


