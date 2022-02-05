function [result_weighted_pca,weightsIm] =  performWeightedPriorReconstruction(dataset,y,param,...
    dim,pilotReconMethods,lambda_for_pilot_recon, baseMethod,lambda_for_base_method,lambda_prior,...
    kk,sliceNumber,outDirectory)

pilotRecons = cell([2,dim]);
for iter = 1:length(pilotReconMethods)
    method = pilotReconMethods(iter);
    switch(method)
        case(0)        
               fname = sprintf('%s/FDK.mat',outDirectory);
               data = load(fname);
               FDK = data.FDK;
               pilotRecons{1} = FDK;
        case(1) 
               lambdaTV = lambda_for_pilot_recon{2};
               fname = sprintf('%s/TV_lambdaTV_%.2f.mat',outDirectory,lambdaTV);
               data = load(fname);
               TV = data.TV;
               pilotRecons{2} = TV;
    end
end

dim = size(FDK);


numCycles = 5;%5;
[result_weighted_pca,weightsIm] =  weighted_prior(dataset,y,pilotReconMethods,lambda_for_pilot_recon, baseMethod,lambda_for_base_method,lambda_prior,dim,param,numCycles,...
                             pilotRecons,kk,sliceNumber,outDirectory);
                   

    
