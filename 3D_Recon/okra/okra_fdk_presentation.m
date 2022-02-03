close all;
clear all;

dataDir = dir('../data/okra/');
numOkras = size(dataDir,1)-2;
outDirectory = 'inputs/presentation';
mkdir(outDirectory);

for i = 1:1;%numOkras
    
    dataDir(i+2).name
    projFileName = sprintf('../data/okra/Preeti_okra%d/projf32_LIN_AP.nc',i);
    ncid = netcdf.open(projFileName);
    proj = double(netcdf.getVar(ncid,1));
    netcdf.close(ncid);
    %------------------------------------------

    startz = 1; endz = size(proj,3); % corresponds to each view
    starty = 100; endy = 255;
    startx = 11; endx = size(proj,1)-10-1; 

    %-------------------------------

    proj = proj(startx:endx, starty: endy, startz:endz); 
    angStepSize = 0.8;


%     mask = find(proj~=-1);
%     minimum = min(proj(mask));
%     proj(mask) = proj(mask) - minimum;
%     proj(proj==-1) = 0;
%     proj = proj./max(proj(:));
%     name = sprintf('okra%d_proj_cropped.avi',i);
%     writevideo (name, proj, 8);

%     proj1 = zeros([168 78 size(proj,3)]);
%     for i = 1:size(proj(:,:,3))
%         img = proj(:,:,i);
%         img = imresize(img,[168 78]);
%         proj1(:,:,i) = img;
%     end
%     proj = proj1;
     size(proj)
    %----------------------------Filtered Back-projection---------------------------

    ParamSetting_okra_pca(angStepSize)
    load('parameters_okra_pca.mat');

    proj_filtered = filtering(proj,param);
    FDK_original = CTbackprojection(proj_filtered, param);
    
    name = sprintf('%s/okra%d_cropped_450views_fdk.mat',outDirectory,i);
    save(name,'FDK_original');
    
    FDK = FDK_original(:,:,1:end-80);
    
    FDK = FDK - min(FDK(:));
    FDK = FDK./max(FDK(:));
    FDKPermuted1 = permute(FDK,[3 1 2]);
    FDKPermuted2 = permute(FDK,[3 2 1]);    
    
    % write to video file
    
    name = sprintf('%s/okra%d_cropped_450views_fdk_z.avi',outDirectory,i);
    writevideo (name, FDK, 4);

    name = sprintf('%s/okra%d_cropped_450views_fdk_y.avi',outDirectory,i);
    writevideo (name, FDKPermuted1, 4);

    name = sprintf('%s/okra%d_cropped_450views_fdk_x.avi',outDirectory,i);
    writevideo (name, FDKPermuted2, 4);
end
















