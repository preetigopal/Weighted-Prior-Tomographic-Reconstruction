function CS = okra_CS(proj,param,dim,lambda0,rel_tol)

Afun = @(z) AConeBeam(z,dim,param);
Atfun = @(z) AtConeBeam(z,dim,param);

y = proj(:) ;

m = size(y,1);
n = param.nx*param.ny*param.nz;       

[x_hat,status,historyNoPrior]=l1_ls_modified(Afun,Atfun,m,n,y,lambda0,rel_tol);
output = reshape(x_hat,[param.nx param.ny param.nz]);
CS = mirt_idctn(output);