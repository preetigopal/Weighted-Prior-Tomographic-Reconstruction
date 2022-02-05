function genEigenSpacePotato_hole(templateNos,outDirectory)

numTemplates = length(templateNos);
numSlices = 100;
sizeOfSlice = [150 150];


volumes  = zeros(sizeOfSlice(1),sizeOfSlice(2),numSlices,numTemplates);

% Read the templates in 3D format-------------------

for i = 1:numTemplates  
   name = sprintf('inputs/potato%d_cropped_900views_fdk.mat',templateNos(i)) ;
   fileData = load(name);
   volumes(:,:,:,i) = fileData.FDK;
   image = volumes(:,:,30,i);
   image = image - min(image(:));
   image = image./max(image(:));
   name = sprintf('inputs/template_%d.png',i);
   imwrite(image,name);
end
 
% store templates and test together------------------
name = sprintf('inputs/potato%d_cropped_900views_fdk.mat',6) ;
fileDataTest = load(name);
testVol = fileDataTest.FDK;
image = testVol(:,:,30);
image = image - min(image(:));
image = image./max(image(:));
name = 'inputs/test.png';
imwrite(image,name);

% fh = figure;imshow(image); hold on;
% rectangle('Position',[55, 78,15 ,15],'EdgeColor','r','LineWidth',2);
% name = sprintf('inputs/test.png');
% f=getframe; imwrite(f.cdata,name);


all_templates = [volumes(:,:,:,1) volumes(:,:,:,2) volumes(:,:,:,3) testVol];
all_templates = all_templates - min(all_templates(:));
all_templates = all_templates./max(all_templates(:));
name = 'templates_and_test';
dataset ='potato';
store_video(all_templates,name,'inputs',dataset);
tempTestImg = [volumes(:,:,30,1) volumes(:,:,30,2) volumes(:,:,30,3) testVol(:,:,30)];
tempTestImg = tempTestImg - min(tempTestImg(:));
tempTestImg = tempTestImg./max(tempTestImg(:));
imwrite(tempTestImg,'inputs/template_test_image.png');
% resize the
% templates--------------------------------------------------------
templates = zeros(sizeOfSlice(1)*sizeOfSlice(2)*numSlices,numTemplates);

for i = 1:numTemplates   
   temp = volumes(:,:,:,i);
   temp = temp(:);   
   templates(:,i) = temp;
end

% Get the templates and compute their mean----------------------------------------------------------

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


% Compute the Covariance Matrix--------------------------------------------------------------------

templates = templates -repmat(meanTemplate,1,numTemplates);
L = templates'*templates;
[W,D] = eig(L);
V = templates*W;
%V = normc(V);
V = normc_fcn(V);
[m n] = size(V);

% picking top k eigen values and their corresponding vectors-----------------------------------------------------
% This forms the eigen space of the covariance matrix of the templates-----------------                  

numDim = numTemplates;%-1;
eigenVals = zeros(1,numDim);
eigenVecs = zeros(m,numDim);

for j = 1:numDim    
    eigenVals(j) = D(n-j+1,n-j+1);
    eigenVecs(:,j) = V(:,n-j+1);
    eigenImage = reshape(eigenVecs(:,j),[sizeOfSlice(1) sizeOfSlice(2) numSlices]);
%     figure;
%     for k = 1:numSlices
%         imshow(eigenImage(:,:,k),[]);
%     end
end
size(eigenVecs)
name = sprintf('%s/EigenSpace_potato_hole.mat',outDirectory);
save(name,'eigenVals','eigenVecs','meanTemplate','minimum','maximum','-v7.3');



