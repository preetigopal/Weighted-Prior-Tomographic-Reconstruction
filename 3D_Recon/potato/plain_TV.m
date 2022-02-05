function result_plain_TV = plain_TV(y,dim,lambdaTV,A_both_functions)

    
addpath('TVReg-master/');

dims=[dim(1) dim(2) dim(3)];
opt_ref.epsb_rel = 1e-4;
opt_ref.k_max    = 20000;
opt_ref.verbose  = 1;
tau = 0.0001;
% Specify nonnegativity constraints
constraint.type = 2;
constraint.c    = 0*ones(prod(dims),1);
constraint.d    = 1*ones(prod(dims),1);
tic;
[x_ref fxk_ref hxk_ref gxk_ref fxkl_ref info_ref] = ...
tvreg_upn(A_both_functions,y,lambdaTV,tau,dims,constraint,opt_ref);
toc;
result_plain_TV = reshape(x_ref,[dim(1) dim(2) dim(3)]);
  
 
    
end
