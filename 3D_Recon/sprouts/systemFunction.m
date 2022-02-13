function output = systemFunction(z,mode,dim,param,m)

if mode==0
    n = dim(1)*dim(2)*dim(3);
    output = [m n];
    
elseif mode==1
    
    X_data = z;
    X = reshape(X_data,[dim(1) dim(2) dim(3)]);
    proj = CTprojection(X,param);
    sizeProj = size(proj);
    b = reshape(proj,[sizeProj(1)*sizeProj(2)*sizeProj(3) 1]);
    output = double(b);
    
elseif mode==2
    
    b_data = z;
    f = load('ProjectionSize.mat');
    b_data = reshape(b_data,[f.sizeProj(1) f.sizeProj(2) f.sizeProj(3)]);
    proj_filtered = filtering(b_data,param);
    backProjImg = CTbackprojection(proj_filtered, param);
    X_final = reshape(backProjImg,[size(backProjImg,1)*size(backProjImg,2)*size(backProjImg,3) 1]);
    output = double(X_final);
end