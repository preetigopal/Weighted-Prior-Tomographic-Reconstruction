function P_pilot = projectPilot_on_EigenSpace_3D(dataset,param,pilotOfTestVol,lambda,A_both_functions,method,sliceNumber,outDirectory)



if strcmp(dataset,'potato6')
    dim = size(pilotOfTestVol);
    templateNos = [1,3,4];
    numTemplates = length(templateNos);
    numSlices = 100;
    sizeOfSlice = [150 150];
    % resize the  templates--------------------------------------------------------
    templates = zeros(sizeOfSlice(1)*sizeOfSlice(2)*numSlices,numTemplates);
    
    % Get the templates and compute their mean----------------------------------------------------------

    for i = 1:numTemplates         
        if strcmp(method,'FDK')       
            fname = sprintf('%s/FDK_template_%d.mat',outDirectory,templateNos(i));
            if(isfile(fname))
                fprintf('Loading FDK reconstructions of template volumes from existing files\n')
                data = load(fname);
                FDK = data.FDK;
            else
                fprintf('Computing FDK reconstructions of template volumes\n');
                projFolderName =  sprintf('inputs/potato%d/projf32_LIN_AP_nc',templateNos(i))
                dirInfo = dir(projFolderName);
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
                %------------------------------------------

                startz = 1; endz = size(proj,3); % corresponds to each view
                starty = 121; endy = size(proj,2)-121+1;
                startx = 125; endx = size(proj,1)-125+1; 

                proj = proj(startx:endx, starty: endy, startz:endz); 
                selectedProj = 1:50:size(proj,3);
                proj = proj(:,:,selectedProj);
                proj1 = zeros(150,150,size(proj,3));
                for ii = 1:size(proj,3)
                    image = imresize(proj(:,:,ii),[150 150]);
                    proj1(:,:,ii) = image;
                end
                proj = proj1;
                FDK = potato_FDK(proj,param,dim);
                save(fname,'FDK');
            end
            temp = FDK - min(FDK(:));
            temp = temp./max(temp(:));
            templates(:,i) = temp(:);

            
        elseif strcmp(method,'TV')   
                fname = sprintf('%s/TV_template_%d.mat',outDirectory,templateNos(i));
            if(isfile(fname))
                fprintf('Loading TV reconstructions of template volumes from existing files\n')
                data = load(fname);
                TV = data.TV;
            else
                fprintf('Computing TV reconstructions of template volumes\n');
                projFolderName =  sprintf('inputs/potato%d/projf32_LIN_AP_nc',templateNos(i));
                dirInfo = dir(projFolderName);
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
                %------------------------------------------
                startz = 1; endz = size(proj,3); % corresponds to each view
                starty = 121; endy = size(proj,2)-121+1;
                startx = 125; endx = size(proj,1)-125+1; 

                proj = proj(startx:endx, starty: endy, startz:endz); 
                selectedProj = 1:50:size(proj,3);
                proj = proj(:,:,selectedProj);
                proj1 = zeros(150,150,size(proj,3));
                for ii = 1:size(proj,3)
                    image = imresize(proj(:,:,ii),[150 150]);
                    proj1(:,:,ii) = image;
                end
                proj = proj1;
                TV = plain_TV(proj(:),dim,lambda,A_both_functions);
                save(fname,'TV');
            end
            temp = TV - min(TV(:));
            temp = temp./max(temp(:));
            templates(:,i) = temp(:);
            
        end  
        
        if i==1
            sumI = zeros(size(templates(:,i)));        
        end
        sumI = sumI + templates(:,i);
    end
    meanTemplate = sumI./numTemplates;


    % Compute the Covariance Matrix--------------------------------------------------------------------

    templates = templates -repmat(meanTemplate,1,numTemplates);
    L = templates'*templates;
    [W,D] = eig(L);
    V = templates*W;
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
       % eigenImage = reshape(eigenVecs(:,j),[sizeOfSlice(1) sizeOfSlice(2) numSlices]);
    %     for k = 1:numSlices
    %         figure;imshow(eigenImage(:,:,k),[]);
    %     end
    end

    %-------------------------------------------------------------------------
    %------------Project the FBP of test image onto this eigen space---------------------------
    %------------------------------------------------------------------------
    % Compute the weights ('alpha') for the FBP of the test image
    % Note: Here- cols of eigenVecs are eigenvectors.

    alpha = eigenVecs'*(pilotOfTestVol(:) - meanTemplate);  


    % Reconstruct the FBP of testIm back from its projection coefficients on the
    % eigen vectors

        coeff = alpha;
        recon = zeros(size(meanTemplate));
        for j = 1:numDim
            recon = recon + (coeff(j)*eigenVecs(:,j));
        end
        P_pilot = recon + meanTemplate;
        P_pilot = reshape(P_pilot,[dim(1) dim(2) dim(3)]);
        %imshow([fdkTestVol(:,:,30) P_fbp(:,:,30)],[]);
        result_projection = [pilotOfTestVol(:,:,sliceNumber) P_pilot(:,:,sliceNumber)];
        result_projection = result_projection - min(result_projection(:));
        result_projection = result_projection./max(result_projection(:));
        if strcmp(method,'FDK') 
            fname = sprintf('%s/projection_on_FDKeigenspace.png',outDirectory);
        elseif strcmp(method,'TV') 
            fname = sprintf('%s/projection_on_TVeigenspace.png',outDirectory);
        end
        imwrite(result_projection,fname);
end