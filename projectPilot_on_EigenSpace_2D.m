function P_pilot = projectPilot_on_EigenSpace_2D(dataset,idx1,pilotOfTestIm,lambda,A_both_functions,method,sliceNo,outDirectory)


% This code is part of the following work which has been submitted to IEEE Transactions on Computational Imaging for peer review.

% Title: "Eliminating object prior-bias from sparse-projection tomographic reconstructions"
% Authors: Preeti Gopal, Sharat Chandran, Imants Svalbe and Ajit Rajwade 

dim = size(pilotOfTestIm);

if strcmp(dataset,'potato')
    templateNos = [1,3,4];
    numTemplates = length(templateNos);
    numSlices = 100;
    sizeOfSlice = [150 150];
    volumes  = zeros(sizeOfSlice(1),sizeOfSlice(2),numSlices,numTemplates);

    for i = 1:numTemplates  
       name = sprintf('../../3D/potato/inputs/potato%d_cropped_900views_fdk.mat',templateNos(i)) ;
       fileData = load(name);
       fdkVol = fileData.FDK;
       resizedVol = fdkVol;%imresize(fdkVol(10:end-13,10:end-13,:),0.5);
       volumes(:,:,:,i) = resizedVol;
    end
end
if strcmp(dataset,'okra')
    templateNos = [4,5,6,7];
    numTemplates = length(templateNos);
    numSlices = 123; % number of slices in a volume;
    sizeOfSlice = [310 310];
    volumes  = zeros(sizeOfSlice(1),sizeOfSlice(2),numSlices,numTemplates);

    for i = 1:numTemplates  
       name = sprintf('../../data/templates/okra/okra6_okra%d_reg_450views_fdk.mat',templateNos(i)) ;
       fileData = load(name);
       fdkVol = fileData.FDK;
       volumes(:,:,:,i) = fdkVol(15:end-14,15:end-14,:);
    end
end

if strcmp(dataset,'sprouts')
    templateNos = [5,6,7,9,10];
    numTemplates = length(templateNos);
    numSlices = 130;%161;
    sizeOfSlice = [260 260];
    volumes  = zeros(sizeOfSlice(1),sizeOfSlice(2),numSlices,numTemplates);

    for i = 1:numTemplates  
       name = sprintf('../../3D/Sprouts/inputs/sprouts%d_givenVolume.mat',templateNos(i)) ;
       fileData = load(name);
       fdkVol = fileData.volume;
       volumes(:,:,:,i) = fdkVol(10:end-4,2:end-2,:);
    end
end

if strcmp(dataset,'tmh_7')
    templateNos = [10,11,12,13,14,15];
    dataFileName = 'RFA2 - Anon/Rfa/RFA 1.5 B30s -';
    numTemplates = length(templateNos);
    numSlices = 43;%161;
    sizeOfSlice = [512 512];
    volumes  = zeros(sizeOfSlice(1),sizeOfSlice(2),numSlices,numTemplates);

    
    for k= 1:numTemplates
        testNo = templateNos(k);
        projFolderName =  sprintf('%s %d/',dataFileName,testNo);
        dirInfo = dir(projFolderName);
        for j = 1:size(dirInfo,1)-2
            fileName = dirInfo(j+2).name;
            filePath = sprintf('%s/%s',projFolderName,fileName);
            image = dicomread(filePath);
            volumes(:,:,j,k) = image;
        end
    end

end

if strcmp(dataset,'tmh_8')
    templateNos = [10,11,12,13,14,15,16];
    dataFileName = 'RFA2 - Anon/Rfa/RFA 1.5 B30s -';
    numTemplates = length(templateNos);
    numSlices = 43;%161;
    sizeOfSlice = [512 512];
    volumes  = zeros(sizeOfSlice(1),sizeOfSlice(2),numSlices,numTemplates);

    
    for k= 1:numTemplates
        testNo = templateNos(k);
        projFolderName =  sprintf('../tmh/%s %d/',dataFileName,testNo);
        dirInfo = dir(projFolderName);
        for j = 1:size(dirInfo,1)-2
            fileName = dirInfo(j+2).name;
            filePath = sprintf('%s/%s',projFolderName,fileName);
            image = dicomread(filePath);
            volumes(:,:,j,k) = image;
        end
    end

end

% Get the templates and compute their mean----------------------------------------------------------
templateIm = cell(1,numTemplates);
for i = 1:numTemplates
    
    in = volumes(:,:,sliceNo,i); 
    if strcmp(dataset,'tmh_7')
        in = in - min(in(:));
        in  = in./max(in(:));
    end
    if strcmp(dataset,'tmh_8')
        in = in - min(in(:));
        in  = in./max(in(:));
    end
    
    temp = in;%double(in);
    radProj = radon(temp,idx1);
    
    if strcmp(method,'FBP')       
        fname = sprintf('%s/FBP_template_%d_num_angles_%d.mat',outDirectory,i,length(idx1));
        if(isfile(fname))
            fprintf('Loading FBP reconstructions of template volumes from existing files')
            data = load(fname);
            fbp = data.fbp;
        else
            fprintf('Computing FBP reconstructions of template volumes');
            fbp = iradon(radProj,idx1,'linear','Cosine');
            fbp = fbp(2:2+dim(1)-1,2:2+dim(2)-1);
            save(fname,'fbp');
        end
        temp = fbp;
    elseif strcmp(method,'TV')   
            fname = sprintf('%s/TV_template_%d_num_angles_%d.mat',outDirectory,i,length(idx1));
        if(isfile(fname))
            fprintf('Loading TV reconstructions of template volumes from existing files')
            data = load(fname);
            TV = data.TV;
        else
            fprintf('Computing TV reconstructions of template volumes');  
            TV = plain_TV(radProj(:),dim,lambda,A_both_functions);
            save(fname,'TV');
        end
        temp = TV;
    end
   
    
 
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
    
    temp = temp - minimum;
    temp = temp./maximum;
    templateIm{i} = temp;
end

meanTemplate = sumI./numTemplates;

% Compute the Covariance Matrix--------------------------------------------------------------------

templates = templates -repmat(meanTemplate,1,numTemplates);
L = templates'*templates;
[W,D] = eig(L);
V = templates*W;
V = normc(V);
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
    %imshow(reshape(eigenVecs(:,j),size(temp)),[]);title(j);pause(0.1);
end


%-------------------------------------------------------------------------
%------------Project the FBP of test image onto this eigen space---------------------------
%------------------------------------------------------------------------
% Compute the weights ('alpha') for the FBP of the test image
% Note: Here- cols of eigenVecs are eigenvectors.
%tic;
alpha = eigenVecs'*(pilotOfTestIm(:) - meanTemplate);  


% Reconstruct the FBP of testIm back from its projection coefficients on the
% eigen vectors

    coeff = alpha;
    recon = zeros(size(meanTemplate));
    for j = 1:numDim
        recon = recon + (coeff(j)*eigenVecs(:,j));
    end
    P_pilot = recon + meanTemplate;
    P_pilot = reshape(P_pilot,[dim(1) dim(2)]);
    %imshow([fbpTestIm P_fbp],[]);
    %toc;
end