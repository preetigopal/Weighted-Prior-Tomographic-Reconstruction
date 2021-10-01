function performWeightedPriorReconstruction(dataset,y,idx1,dim,pilotReconMethods,lambda_for_pilot_recon, baseMethod,lambda_for_base_method,lambda_prior_list,kk,sliceNum,outDirectory)

pilotRecons = cell([2,dim]);
for iter = 1:length(pilotReconMethods)
    method = pilotReconMethods(iter);
    switch(method)
        case(0)        
               fname = sprintf('%s/fbp.mat',outDirectory);
               data = load(fname);
               result_fbp = data.result_fbp;
               pilotRecons{1} = result_fbp;
        case(1) 
               lambdaTV = lambda_for_base_method;
               fname = sprintf('%s/tv_lambdaTV_%.3f.mat',outDirectory,lambdaTV);
               data = load(fname);
               result_tv = data.result_tv;
               pilotRecons{2} = result_tv;
    end
end

dim = size(result_fbp);

for iter = 1:length(lambda_prior_list)
    lambda_prior = lambda_prior_list(iter);
    numCycles = 5;
    [result_weighted_pca,weightsIm] =  weighted_prior(dataset,y,pilotReconMethods,lambda_for_pilot_recon, baseMethod,lambda_for_base_method,lambda_prior,dim,idx1,numCycles,...
                                 pilotRecons,sliceNum,outDirectory,kk);

    fname = sprintf('%s/weighted_prior_kk_%d_lambda_prior_%f.mat',outDirectory,kk,lambda_prior);
    save(fname,'result_weighted_pca');
    name = sprintf('%s/weighted_prior_kk_%d_lambda_prior_%f.png',outDirectory,kk,lambda_prior);
    temp = result_weighted_pca;
    temp = temp - min(temp(:));
    temp = temp/max(temp(:));
    imwrite(temp,name);    

    fname = sprintf('%s/weightsIm_kk_%d_lambda_prior_%f.mat',outDirectory,kk,lambda_prior);
    save(fname,'weightsIm');
    name = sprintf('%s/weightsIm_kk_%d_lambda_prior_%f.png',outDirectory,kk,lambda_prior);
    temp = weightsIm;
    temp = temp - min(temp(:));
    temp = temp/max(temp(:));
    imwrite(temp,name);  
end
