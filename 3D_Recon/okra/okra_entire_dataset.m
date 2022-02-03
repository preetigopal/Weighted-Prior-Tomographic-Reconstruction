close all;
clear all;

templateNos = [1,2,3,4,5,6, 7];
numTemplates = length(templateNos);
numSlices = 123; % number of slices in a volume;
sizeOfSlice = [338 338];


volumes  = zeros(sizeOfSlice(1),sizeOfSlice(2),numSlices,numTemplates);

% Read the templates in 3D format-------------------

for i = 1:numTemplates  
   name = sprintf('inputs/okra%d_okra6_reg_cropped_450views_fdk.mat',templateNos(i)) ;
   fileData = load(name);
   volumes(:,:,:,i) = fileData.FDK;
   
    temp = volumes(:,:,30,i);
    fname = sprintf('slice%d.mat',i);
    save(fname,'temp');    
    temp = temp - min(temp(:));
    temp = temp./max(temp(:));    
    name = sprintf('slice%d.png',i);
    imwrite(temp,name);
end
 
allVolumes = [volumes(:,:,1:60,1) volumes(:,:,1:60,2) volumes(:,:,1:60,3) volumes(:,:,1:60,4) volumes(:,:,1:60,5) volumes(:,:,1:60,6) volumes(:,:,1:60,7)];
allVolumes = allVolumes - min(allVolumes(:));
allVolumes = allVolumes./max(allVolumes(:));
name = 'entire_dataset';
dataset ='okra';
store_video(allVolumes,name,'inputs',dataset);
tempTestImg = [volumes(:,:,30,1) volumes(:,:,30,2) volumes(:,:,30,3) volumes(:,:,30,4) volumes(:,:,30,5) volumes(:,:,30,6) volumes(:,:,30,7) ];
tempTestImg = tempTestImg - min(tempTestImg(:));
tempTestImg = tempTestImg./max(tempTestImg(:));
imwrite(tempTestImg,'inputs/okra_entire_dataset.png');


