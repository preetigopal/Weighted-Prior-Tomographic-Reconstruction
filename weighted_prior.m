function [result_weighted_prior,W] =  weighted_prior(dataset,y_given,pilotReconMethods,lambda_for_pilot_recon, baseMethod,lambda_for_base_method,lambda_prior,dim,idx1,numCycles,...
                                 pilotRecons,sliceNum,outDirectory,kk)
                             

% This code is part of the following work which has been submitted to Transactions of Computational Imaging for peer review.

% Title: "Eliminating object prior-bias from sparse-projection tomographic reconstructions"
% Authors: Preeti Gopal, Sharat Chandran, Imants Svalbe and Ajit Rajwade 

% This function computes the weighted prior reconstruction by first
% estimating the weights-map W and then computing the reconstruction by
% alternate minimization (of Eq.4 in the paper).


if strcmp(dataset,'tmh_7')
    templateNos = [10,11,12,13,14,15];
    dataFileName = 'RFA2 - Anon/Rfa/RFA 1.5 B30s -';
    numSlices = 43;%161;
    sizeOfSlice = [512 512];
%     sliceNum = 27;
    
    [eigenVecs,meanTemplateOriginal] = genEigenSpace2Dtmh(templateNos,sliceNum,dataFileName,numSlices,sizeOfSlice,outDirectory);

    
elseif strcmp(dataset,'tmh_8')
    templateNos = [10,11,12,13,14,15,16]
    dataFileName = 'RFA2 - Anon/Rfa/RFA 1.5 B30s -';
    numSlices = 43;%161;
    sizeOfSlice = [512 512];
%     sliceNum = 27;

    [eigenVecs,meanTemplateOriginal] = genEigenSpace2Dtmh(templateNos,sliceNum,dataFileName,numSlices,sizeOfSlice,outDirectory);
else
    fname = sprintf('EigenSpace_2D_%s.mat',dataset); 
    es = load(fname);
    eigenVecs = es.eigenVecs;
    meanTemplateOriginal = es.meanTemplate;
end
numAngles = length(idx1);
m = size(y_given,1);
A_both_functions =@(z,mode) systemFunction(z,mode,idx1,dim,numAngles,m);    

%%  --------------------------construct the weight map---------------------------------------

for iter = 1:length(pilotReconMethods)
    method = pilotReconMethods(iter);
    switch(method)
        case(0)
            lambda = 1; % dummy value
            FBP = pilotRecons{1};
            method_name = 'FBP';
            P_FBP = projectPilot_on_EigenSpace_2D(dataset,idx1,FBP,lambda,A_both_functions,method_name,sliceNum,outDirectory);
            fname = sprintf('%s/ProjectedTestFBPIm.mat',outDirectory);
            save(fname,'P_FBP');
            name = sprintf('%s/FBP_P_FBP.png',outDirectory);
            temp = [FBP P_FBP];
            %figure;imshow(temp,[]);
            temp = temp - min(temp(:));
            temp = temp/max(temp(:));
            imwrite(temp,name);
        case(1)
            lambda_for_method = lambda_for_pilot_recon(2); 
            lambda_TV = lambda_for_method{1};
            TV = pilotRecons{2};
            method_name = 'TV';
            P_TV = projectPilot_on_EigenSpace_2D(dataset,idx1,TV,lambda_TV,A_both_functions,method_name,sliceNum,outDirectory);
            fname = sprintf('%s/ProjectedTestTVIm.mat',outDirectory);
            save(fname,'P_TV');
            name = sprintf('%s/TV_P_TV.png',outDirectory);
            temp = [TV P_TV];
            %figure;imshow(temp,[]);
            temp = temp - min(temp(:));
            temp = temp/max(temp(:));
            imwrite(temp,name);  
    end
end

k_value = kk;
%tic;
error_fbp = abs(P_FBP - FBP);
error_tv = abs(P_TV - TV);
error_all = bsxfun(@min,error_fbp,error_tv);
W = (1./(1+k_value.*error_all));
%toc;

%% ---------------begin reconstruction----------------------------------
addpath('TVReg-master/');
meanTemplate = meanTemplateOriginal;
output = zeros(dim);
alphas = cell(1,numCycles);    % assume initial set of alphas to be all zeros
tic;
for numIter = 1:numCycles  

    display (numIter)

    if numIter==1   
        
        y = y_given;        
        priorIm = meanTemplate.*W(:);
        temp = lambda_prior.*priorIm;  % including the prior term
        y = cat(1,y,temp);

    else
          
        y = y_given;

        result = output;
        result  = result(:);
        % alphas{numIter} =  eigenVecs'*(result- meanTemplate);
        % Start: compute alphas---------------
        weightedEigenVecs = eigenVecs.*W(:);
        matrix = (weightedEigenVecs)' * (weightedEigenVecs);
        z = W(:).*(result - meanTemplate);
        alphas{numIter} = matrix\(eigenVecs'*(W(:).*z));
        
        % End: compute alphas------------------
        
        priorIm = (meanTemplate + eigenVecs*alphas{numIter}).*W(:);        

        temp = lambda_prior.*priorIm;
        y = cat(1,y,temp);

    end 
    
    if strcmp(baseMethod,'TV')
        m = size(y,1);
        numAngles = length(idx1);
        A_both_functions_prior =@(z,mode) systemFunction_prior(z,mode,idx1,dim,numAngles,m,lambda_prior,W);
        
        dims=[dim(1) dim(2)];
        opt_ref.epsb_rel = 1e-4;
        opt_ref.k_max    = 40000;
        opt_ref.verbose  = 1;
        tau = 0.0001;
        % Specify nonnegativity constraints
        constraint.type = 2;
        constraint.c    = 0*ones(prod(dims),1);
        constraint.d    = 1*ones(prod(dims),1);
        [x_ref fxk_ref hxk_ref gxk_ref fxkl_ref info_ref] = ...
        tvreg_upn(A_both_functions_prior,y,lambda_for_base_method,tau,dims,constraint,opt_ref);
        output = reshape(x_ref,[dim(1) dim(2)]);
    end
  

end
toc;
result_weighted_prior = output;    