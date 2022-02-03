close all;
clear all;
dataset = 'okra3';

projFileName = sprintf('../data/okra/Preeti_%s/projf32_LIN_AP.nc',dataset);
ncid = netcdf.open(projFileName);
proj = double(netcdf.getVar(ncid,1));
netcdf.close(ncid);
%------------------------------------------
  
startz = 1; endz = size(proj,3); % corresponds to each view
starty = 100; endy = 255;%size(proj,2)-1;
startx = 11; endx = size(proj,1)-10-1; 

proj = proj(startx:endx, starty: endy, startz:endz);
selectedProj = 1:10:size(proj,3);
proj = proj(:,:,selectedProj);
angStepSize = 8;
ParamSetting_okra_pca(angStepSize)
load('parameters_okra_pca.mat');

dim = [param.nx param.ny param.nz]; 
%----------------------------Filtered Back-projection---------------------------

proj_filtered = filtering(proj,param);
FDK = CTbackprojection(proj_filtered, param);

FDK = FDK - min(FDK(:));
FDK= FDK./max(FDK(:));
FDKPermuted1 = permute(FDK,[3 1 2]);
FDKPermuted2 = permute(FDK,[3 2 1]);

name = sprintf('%s_reg_45views_fdk.mat',dataset);
% FDKfile = load(name);
% FDK = FDKfile.FDK;
save(name,'FDK');
% % 
% % % write to video file
% % 
name = sprintf('%s_reg_45views_fdk_z.avi',dataset);
writevideo (name, FDK,8);
name = sprintf('%s_reg_45views_fdk_y.avi',dataset);
writevideo (name, FDKPermuted1, 8);
name = sprintf('%s_reg_45views_fdk_x.avi',dataset);
writevideo (name, FDKPermuted2, 8);
% 
% %-------------------------Simple CS (No prior)-------------------------------
% 
lambda0 = 1;
Afun = @(z) AConeBeam(z,dim,param);
Atfun = @(z) AtConeBeam(z,dim,param);

y = proj(:) ;
% 
rel_tol = 0.1; % relative target duality gap
m = size(y,1);
n = param.nx*param.ny*param.nz;       
%tic;
[x_hat,status,historyNoPrior]=l1_ls_modified(Afun,Atfun,m,n,y,lambda0,rel_tol);
output = reshape(x_hat,[param.nx param.ny param.nz]);
CS = mirt_idctn(output); 
%toc;

CS = CS - min(CS(:));
CS = CS./max(CS(:));
CSPermuted1 = permute(CS,[3 1 2]);
CSPermuted2 = permute(CS,[3 2 1]);

name = sprintf('%s_reg_45views_CS.mat',dataset);
% CSfile = load(name);
% CS = CSfile.CS;
save(name,'CS');

name = sprintf('%s_reg_45views_CS_z.avi',dataset);
writevideo (name, CS, 8);
name = sprintf('%s_reg_45views_CS_y.avi',dataset);
writevideo (name, CSPermuted1, 8);
name = sprintf('%s_reg_45views_CS_x.avi',dataset);
writevideo (name, CSPermuted2, 8);
% % 
%-------------------------------Reconstructing with prior--------------------------------------------------------------------
lambda1_values = 0.25;
numCycles = 3;
totalNumIterations = length(lambda1_values);
primaryObj = zeros(1,totalNumIterations);  
sizeOfLoopIndxMatrix = length(lambda1_values);

es = load('EigenSpace_okra_hole.mat');
eigenVecs = es.eigenVecs;
eigenVals = es.eigenVals;
meanTemplate = es.meanTemplate;  

name = sprintf('%s_cropped_450views_fdk.mat',dataset);
testVolFile = load(name);
testVol = testVolFile.FDK;

% resizedMeanTemplate = reshape(meanTemplate,size(testVol));
% % meanTempProj = CTprojection(resizedMeanTemplate,param);
% % meanTempProj_filtered = filtering(meanTempProj,param);
% % FDK_meanTemplate = CTbackprojection(meanTempProj_filtered, param);
% % 
% % FDK_meanTemplate = FDK_meanTemplate - min(FDK_meanTemplate(:));
% % FDK_meanTemplate = FDK_meanTemplate./max(FDK_meanTemplate(:));
% % save('FDK_meanTemplate.mat','FDK_meanTemplate');
% FDKmeanFile = load('FDK_meanTemplate.mat');
% FDK_meanTemplate = FDKmeanFile.FDK_meanTemplate;

% meanTempProj = CTprojection(resizedMeanTemplate,param);
% meanTempProj_filtered = filtering(meanTempProj,param);
% FDK_meanTemplate = CTbackprojection(meanTempProj_filtered, param);
% 
% FDK_meanTemplate = FDK_meanTemplate - min(FDK_meanTemplate(:));
% FDK_meanTemplate = FDK_meanTemplate./max(FDK_meanTemplate(:));
% save('FDK_meanTemplate.mat','FDK_meanTemplate');
% FDKmeanFile = load('FDK_meanTemplate.mat');
% FDK_meanTemplate = FDKmeanFile.FDK_meanTemplate;


% 
% tic; 
for ix = 1:totalNumIterations

    [lamb] = ind2sub(sizeOfLoopIndxMatrix,ix);        
    lambda1 = lambda1_values(lamb);        
    output = zeros(size(testVol));

    %------------------------------------------reconstruction with l2 norm of prior-------------------------------------------

    alphas = cell(1,numCycles);    % assume initial set of alphas to be all zeros

                for numIter = 1:numCycles
                    
                    display (numIter)

                    if numIter==1

                        y = proj(:) ;

                        priorVol = meanTemplate;
                        priorVol = reshape(priorVol,size(testVol));
                        
                        templateVolFile = load('okra6_cropped_45views_fdk');
                        templateVol = templateVolFile.FDK;
                        [optimizer,metric] = imregconfig('monomodal');
%                         file = load('okra4_okra6_tform.mat');
%                         InitialTform = file.tform;
                        tform = imregtform(templateVol(38:308,138:218,1:55),FDK(38:308,138:218,1:55),'rigid',optimizer,metric);%'InitialTransformation',InitialTform);
                        priorVol = imwarp(priorVol,tform,'OutputView',imref3d(size(FDK)));                        
                        
                        priorVol = priorVol(:);

                        temp = lambda1*priorVol;  % including the prior term
                        y = cat(1,y,temp);


                    else

                        y = proj(:) ;

                        result = output;
                        result  = result(:);
                        alphas{numIter} =  eigenVecs'*(result- meanTemplate);
                        priorVol = meanTemplate + eigenVecs*alphas{numIter};
                        priorVol = reshape(priorVol,size(testVol));
                        
 
                        %tform = imregtform(priorVol(:,:,1:55),output(:,:,1:55),'rigid',optimizer,metric);
                        priorVol = imwarp(priorVol,tform,'OutputView',imref3d(size(output)));             
                        
                        priorVol = priorVol(:);

                        temp = lambda1*priorVol;
                        y = cat(1,y,temp);

                    end                         

                    m = size(y,1);
                    n = param.nx*param.ny*param.nz;
                    Afun = @(z) ARadon1_3D(z,dim,lambda1,param);
                    Atfun = @(z) AtRadon1_3D(z,dim,lambda1,param);

                    [x_hat,status,history]=l1_ls_modified(Afun,Atfun,m,n,y,lambda0,rel_tol);

                    startPt = 1;
                    endPt = size(x_hat,1);
                    xp = x_hat(startPt:endPt);

                    output = reshape(xp, [dim(1) dim(2) dim(3)]);
                    output = mirt_idctn(output); 
                    
                    priorVol = reshape(priorVol,size(testVol));
                    priorVol = priorVol - min(priorVol(:));
                    priorVol = priorVol./max(priorVol(:));
                    priorVolPermuted1 = permute(priorVol,[3 1 2]);
                    priorVolPermuted2 = permute(priorVol,[3 2 1]);
                    name = sprintf('priorVol%d_z.avi',numIter);
                    writevideo (name,priorVol,8);
                    name = sprintf('priorVol%d_y.avi',numIter);
                    writevideo (name,priorVolPermuted1,8);
                    name = sprintf('priorVol%d_x.avi',numIter);
                    writevideo (name,priorVolPermuted2,8);
                    
                    name = sprintf('tform%d.mat',numIter);
                    save(name,'tform');
                    
%                     FDK_meanTemplate = reshape(FDK_meanTemplate,size(testVol));
%                     FDK_meanTemplate = FDK_meanTemplate - min(FDK_meanTemplate(:));
%                     FDK_meanTemplate = FDK_meanTemplate./max(FDK_meanTemplate(:));
%                     FDK_meanTemplatePermuted1 = permute(FDK_meanTemplate,[3 1 2]);
%                     FDK_meanTemplatePermuted2 = permute(FDK_meanTemplate,[3 2 1]);
%                     name = sprintf('FDK_meanTemplate%d_z.avi',numIter);
%                     writevideo (name,FDK_meanTemplate,8);
%                     name = sprintf('FDK_meanTemplate%d_y.avi',numIter);
%                     writevideo (name,FDK_meanTemplatePermuted1,8);
%                     name = sprintf('FDK_meanTemplate%d_x.avi',numIter);
%                     writevideo (name,FDK_meanTemplatePermuted2,8);
                    
                    intermediateOutput = reshape(output,size(testVol));
                    intermediateOutput = intermediateOutput - min(intermediateOutput(:));
                    intermediateOutput = intermediateOutput./max(intermediateOutput(:));
                    intermediateOutputPermuted1 = permute(intermediateOutput,[3 1 2]);
                    intermediateOutputPermuted2 = permute(intermediateOutput,[3 2 1]);
                    name = sprintf('intermediateOutput%d_z.avi',numIter);
                    writevideo (name,intermediateOutput,8);
                    name = sprintf('intermediateOutput%d_y.avi',numIter);
                    writevideo (name,intermediateOutputPermuted1,8);
                    name = sprintf('intermediateOutput%d_x.avi',numIter);
                    writevideo (name,intermediateOutputPermuted2,8);
                    
                end
    primaryObj(ix) = history(2,end);      
    pca = output;    
% 
%---------------------compare -----------------------------      
end
%toc;

pca = pca - min(pca(:));
pca = pca./max(pca(:));
pcaPermuted1 = permute(pca,[3 1 2]);
pcaPermuted2 = permute(pca,[3 2 1]);
name = sprintf('%s_reg_45views_pca.mat',dataset);
save(name,'pca');

name = sprintf('%s_reg_45views_pca_z.avi',dataset);
writevideo (name, pca, 8);
name = sprintf('%s_reg_45views_pca_y.avi',dataset);
writevideo (name, pcaPermuted1, 8);
name = sprintf('%s_reg_45views_pca_x.avi',dataset);
writevideo (name, pcaPermuted2, 8);

%%-----------------Computing error metric--------------------------
mse = mean((FDK - testVol).^2,3);
mse = mean2(mse);
msi = mean(testVol.^2,3);
msi = mean2(msi);
relMseValFDK = mse/msi;
ssimFDK = ssim(double(FDK),double(testVol));

mse = mean((CS - testVol).^2,3);
mse = mean2(mse);
msi = mean(testVol.^2,3);
msi = mean2(msi);
relMseVal_plainCS = mse/msi;
ssim_plainCS = ssim(double(CS),double(testVol));

mse = mean((pca - testVol).^2,3);
mse = mean2(mse);
msi = mean(testVol.^2,3);
msi = mean2(msi);
relMseVal_global_prior = mse/msi;
ssim_globalPrior = ssim(double(pca),double(testVol));

save('resultOkra_hole.mat','testVol','FDK','CS','pca','relMseValFDK','relMseVal_plainCS','relMseVal_global_prior','ssimFDK','ssim_plainCS','ssim_globalPrior');







