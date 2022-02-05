function [result_weighted_prior,W] =  weighted_prior(dataset,y_given,pilotReconMethods,lambda_for_pilot_recon, baseMethod,lambda_for_base_method,lambda_prior,dim,param,numCycles,...
                                 pilotRecons,kk,sliceNumber,outDirectory)

fname = sprintf('%s/EigenSpace_%s_hole.mat',outDirectory,'potato'); 
es = load(fname);
eigenVecs = es.eigenVecs;
meanTemplateOriginal = es.meanTemplate;
m = size(y_given,1);
A_both_functions =@(z,mode) systemFunction(z,mode,dim,param,m);  

%%  --------------------------construct the weight map---------------------------------------

for iter = 1:length(pilotReconMethods)
    method = pilotReconMethods(iter);
    switch(method)
        case(0)
            lambda = 1; % dummy value
            FDK = pilotRecons{1};
            FDK = FDK - min(FDK(:));
            FDK = FDK./max(FDK(:));	
            method_name = 'FDK';
            P_FDK = projectPilot_on_EigenSpace_3D(dataset,param,FDK,lambda,A_both_functions,method_name,sliceNumber,outDirectory);
            fname = sprintf('%s/ProjectedTestFDK.mat',outDirectory);
            save(fname,'P_FDK');
            name = sprintf('%s/FDK_P_FDK.png',outDirectory);
            temp = [FDK(:,:,sliceNumber) P_FDK(:,:,sliceNumber)];
            %figure;imshow(temp,[]);
            temp = temp - min(temp(:));
            temp = temp/max(temp(:));
            imwrite(temp,name);
        case(1)
            lambda_for_method = lambda_for_pilot_recon(2); 
            lambda_TV = lambda_for_method{1};
            TV = pilotRecons{2};
            TV = TV - min(TV(:));
            TV = TV./max(TV(:));
            method_name = 'TV';
            P_TV = projectPilot_on_EigenSpace_3D(dataset,param,TV,lambda_TV,A_both_functions,method_name,sliceNumber,outDirectory);
            fname = sprintf('%s/ProjectedTestTV.mat',outDirectory);
            save(fname,'P_TV');
            name = sprintf('%s/TV_P_TV.png',outDirectory);
            temp = [TV(:,:,sliceNumber) P_TV(:,:,sliceNumber)];
            %figure;imshow(temp,[]);
            temp = temp - min(temp(:));
            temp = temp/max(temp(:));
            imwrite(temp,name);  
    end
end

k_value = kk;
error_fbp = abs(P_FDK - FDK);
error_tv = abs(P_TV - TV);
error_all = bsxfun(@min,error_fbp,error_tv);
W = (1./(1+k_value.*error_all));
W = double(W);
% W = W - min(W(:));
% W = W./max(W(:));

meanTemplateReshaped = reshape(meanTemplateOriginal,dim);

fname = sprintf('%s/meanTemplateReshaped',outDirectory);
save(fname,'meanTemplateReshaped');
name = sprintf('%s/meanTemplateReshaped.png',outDirectory);
temp = meanTemplateReshaped(:,:,sliceNumber);
temp = temp - min(temp(:));
temp = temp./max(temp(:));
imwrite(temp,name);

meanTemplateReshaped1 = meanTemplateReshaped - min(meanTemplateReshaped(:));
meanTemplateReshaped1 = meanTemplateReshaped1./max(meanTemplateReshaped1(:));
summary_video_name = sprintf('meanTemplate'); 
store_video(meanTemplateReshaped1,summary_video_name,outDirectory,'potato'); 

%% ---------------begin reconstruction----------------------------------

meanTemplate = meanTemplateOriginal;
output = zeros(dim);
alphas = cell(1,numCycles);    % assume initial set of alphas to be all zeros

for numIter = 1:numCycles  

    display (numIter)

    if numIter==1   
        
        y = y_given;        
        priorIm = meanTemplate;
        temp = lambda_prior.*priorIm.*W(:);  % including the prior term
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
        
        priorIm = (meanTemplate + eigenVecs*alphas{numIter});        

        temp = lambda_prior.*priorIm.*W(:);
        y = cat(1,y,temp);

    end 
    
    if strcmp(baseMethod,'TV')
        m = size(y,1);
        A_both_functions_prior =@(z,mode) systemFunction_prior(z,mode,param,dim,m,lambda_prior,W);
        addpath('TVReg-master/');
        opt_ref.epsb_rel = 1e-4;
        opt_ref.k_max    = 200;
        opt_ref.verbose  = 1;
        tau = 0.0001;
        % Specify nonnegativity constraints
        constraint.type = 2;
        constraint.c    = 0*ones(prod(dim),1);
        constraint.d    = 1*ones(prod(dim),1);
        [x_ref fxk_ref hxk_ref gxk_ref fxkl_ref info_ref] = ...
        tvreg_upn(A_both_functions_prior,y,lambda_for_base_method,tau,dim,constraint,opt_ref);
        output = reshape(x_ref,[dim(1) dim(2) dim(3)]);
    end
    priorImReshaped = reshape(priorIm,dim);
    fname = sprintf('%s/priorIm_iteration_%d_lambda_prior_%f.png',outDirectory,numIter,lambda_prior);
    save(fname,'priorImReshaped');
    name = sprintf('%s/priorIm_iteration_%d_lambda_prior_%f.png',outDirectory,numIter,lambda_prior);
    temp = priorImReshaped(:,:,sliceNumber);
    temp = temp - min(temp(:));
    temp = temp./max(temp(:));
    imwrite(temp,name);
    
    priorImReshaped1 = priorImReshaped - min(priorImReshaped(:));
    priorImReshaped1 = priorImReshaped1./max(priorImReshaped1(:));
    summary_video_name = sprintf('priorIm_iteration_%d_lambda_prior_%f',numIter,lambda_prior); 
    store_video(priorImReshaped1,summary_video_name,outDirectory,'potato');  
    
    fname = sprintf('%s/weighted_prior_output_iteration_%d_lambda_prior_%f.png',outDirectory,numIter,lambda_prior);
    save(fname,'output');
    name = sprintf('%s/weighted_prior_output_iteration_%d_lambda_prior_%f.png',outDirectory,numIter,lambda_prior);
    temp = output(:,:,sliceNumber);
    temp = temp - min(temp(:));
    temp = temp./max(temp(:));
    imwrite(temp,name);
    
    output1 = output - min(output(:));
    output1 = output1./max(output1(:));
    summary_video_name = sprintf('wighted_prior_output_iteration_%d_lambda_prior_%f',numIter,lambda_prior); 
    store_video(output1,summary_video_name,outDirectory,'potato');  


end
  
result_weighted_prior = output;  