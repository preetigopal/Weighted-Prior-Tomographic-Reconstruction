function output = systemFunction(z,mode,idx,dim,numAngles,m)

% This code is part of the following work which has been submitted to Transactions of Computational Imaging for peer review.

% Title: "Eliminating object prior-bias from sparse-projection tomographic reconstructions"
% Authors: Preeti Gopal, Sharat Chandran, Imants Svalbe and Ajit Rajwade 

if mode==0
    n = dim(1)*dim(2);
    output = [m n];
    
elseif mode==1
    
    X_data = z;
    startPt = 1;
    endPt = dim(1)*dim(2);
    X = X_data(startPt:endPt);
    X = reshape(X,[dim(1) dim(2)]);

    startPt = 1;
    endPt = numAngles;
    angles = idx(startPt:endPt);

    radProj = radon(X,angles);  

    b1 = reshape(radProj,[size(radProj,1)*size(radProj,2) 1]);

    b = b1;
    output = b;
    
elseif mode==2
    
    b_data = z;
    s = size(b_data,1);

    startPt = 1;
    endPt = s;
    b = b_data(startPt:endPt);
    b = reshape(b,[(size(b,1))/numAngles numAngles]);

    startPt = 1;
    endPt = numAngles;
    angles = idx(startPt:endPt);    

    backProjImg  =  iradon(b,angles,'linear','Cosine');

    backProjImg = backProjImg(2:2+dim(1)-1,2:2+dim(2)-1);
    X = reshape(backProjImg,[size(backProjImg,1)*size(backProjImg,2) 1]);       


    X_final = X;
    output = X_final;
end