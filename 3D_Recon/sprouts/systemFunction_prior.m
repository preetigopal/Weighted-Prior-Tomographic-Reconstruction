function output = systemFunction_prior(z,mode,param,dim,m,lambda_prior,W)

if mode==0
    n = dim(1)*dim(2)*dim(3);
    output = [m n];
    
elseif mode==1  %forward system matrix operator
    
    X = z;
    X = reshape(X,[dim(1) dim(2) dim(3)]);
    
    proj = CTprojection(X,param);
    sizeProj = size(proj);
    b1 = reshape(proj,[sizeProj(1)*sizeProj(2)*sizeProj(3) 1]);
    save('ProjectionSize.mat','sizeProj');
    
    vectX =  lambda_prior.*reshape(X,[dim(1)*dim(2)*dim(3) 1]).*W(:); 
    output = cat(1,b1,vectX);
    output = double(output);
    
elseif mode==2
    
    b_data = z;
    s1 = size(b_data,1)-(dim(1)*dim(2)*dim(3)); %m
    s2 = size(b_data,1);  %(m+n)

    startPt = 1;
    endPt =  s1;
    b = b_data(startPt:endPt);
    
    startPt = endPt + 1;
    endPt = s2;
    Y2 = b_data(startPt:endPt);

    f = load('ProjectionSize.mat');
    b = reshape(b,[f.sizeProj(1) f.sizeProj(2) f.sizeProj(3)]);
    proj_filtered = filtering(b,param);
    backProjImg = CTbackprojection(proj_filtered, param);

    X = reshape(backProjImg,[size(backProjImg,1)*size(backProjImg,2)*size(backProjImg,3) 1]); 
       
%-------------- Modification for the PICCS part -starts------------------
    X = X + (lambda_prior).*Y2.*W(:);

%-------------- Modification for the PICCS part -ends------------------

    output = double(X);

end