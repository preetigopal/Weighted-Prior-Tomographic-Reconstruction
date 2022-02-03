function P_pilot = projectPilot_on_EigenSpace_3D(dataset,param,pilotOfTestVol,lambda,A_both_functions,method,sliceNumber,outDirectory)


if strcmp(dataset,'okra3')
    dim = size(pilotOfTestVol);
    templateNos = [4,5,6,7];
    numTemplates = length(templateNos);
    numSlices = 123; 
    sizeOfSlice = [338 338];
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
                name = sprintf('inputs/okra%d_okra6_reg_cropped_450views_fdk.mat',templateNos(i))
                fileData  = load(name);
                FDK_high_quality = fileData.FDK;
                proj = double(CTprojection(FDK_high_quality,param));
                
                % START: If raw projections of templates are to be used----------
                
%                 projFileName = sprintf('inputs/Preeti_okra%d/projf32_LIN_AP.nc',startTemplateNum + (i-1));
%                 ncid = netcdf.open(projFileName);
%                 proj = double(netcdf.getVar(ncid,1));
%                 netcdf.close(ncid);
% 
%                 startz = 1; endz = size(proj,3); % corresponds to each view
%                 starty = 100; endy = 255;
%                 startx = 11; endx = size(proj,1)-10-1; 
% 
%                 proj = proj(startx:endx, starty: endy, startz:endz);
%                 selectedProj = 1:10:size(proj,3);
%                 proj = proj(:,:,selectedProj);

                % END: If raw projections of templates are to be used----------

                FDK = okra_FDK(proj,param,dim);
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
                name = sprintf('inputs/okra%d_okra6_reg_cropped_450views_fdk.mat',templateNos(i))
                fileData = load(name);
                FDK_high_quality = fileData.FDK;
                proj = double(CTprojection(FDK_high_quality,param));
                
                % START: If raw projections of templates are to be used----------
                
%                 projFileName = sprintf('inputs/Preeti_okra%d/projf32_LIN_AP.nc',startTemplateNum + (i-1));
%                 ncid = netcdf.open(projFileName);
%                 proj = double(netcdf.getVar(ncid,1));
%                 netcdf.close(ncid);
%             
%                 startz = 1; endz = size(proj,3); % corresponds to each view
%                 starty = 100; endy = 255;
%                 startx = 11; endx = size(proj,1)-10-1; 
%             
%                 proj = proj(startx:endx, starty: endy, startz:endz);
%                 selectedProj = 1:10:size(proj,3);
%                 proj = proj(:,:,selectedProj);

                % END: If raw projections of templates are to be used----------

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