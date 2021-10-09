function performPilotReconstruction(y,idx1,dim,pilotReconMethods,lambda_for_pilot_recon,outDirectory)

% This code is part of the following work which has been submitted to IEEE Transactions on Computational Imaging for peer review.

% Title: "Eliminating object prior-bias from sparse-projection tomographic reconstructions"
% Authors: Preeti Gopal, Sharat Chandran, Imants Svalbe and Ajit Rajwade 


% This function computes FBP and TV reconstructions of the 2D data.

numAngles = length(idx1)
m = size(y,1);
A_both_functions =@(z,mode) systemFunction(z,mode,idx1,dim,numAngles,m);

for iter = 1:length(pilotReconMethods)
    method = pilotReconMethods(iter);
    lambda_for_the_method = lambda_for_pilot_recon{iter};
    
    switch(method)
        case(0)
            fname = sprintf('%s/fbp.mat',outDirectory);
            if(isfile(fname))
                continue;
            else
                %tic;
                result_fbp = filteredBackProjection2D(y,numAngles,idx1,dim);
                %toc;
                save(fname,'result_fbp');
                name = sprintf('%s/fbp.png',outDirectory);
                temp = result_fbp;
                temp = temp - min(temp(:));
                temp = temp/max(temp(:));
                imwrite(temp,name);
            end
        
        case(1)            
            lambdaTV_list =lambda_for_the_method;
            for i = 1:length(lambdaTV_list)
                lambdaTV = lambdaTV_list(i)
                %fname = sprintf('%s/tv_lambdaTV_%.2f.mat',outDirectory,lambdaTV);
                fname = sprintf('%s/tv_lambdaTV_%.3f.mat',outDirectory,lambdaTV); 
                if(isfile(fname))
                    continue;
                else
                    %tic;
                    result_tv = plain_TV(y,dim,lambdaTV,A_both_functions);
                    %toc;
                    save(fname,'result_tv');
                    name = sprintf('%s/tv_lambdaTV_%.3f.png',outDirectory,lambdaTV);
                    temp = result_tv;
                    temp = temp - min(temp(:));
                    temp = temp/max(temp(:));
                    imwrite(temp,name);
                end
            end
    end
end

