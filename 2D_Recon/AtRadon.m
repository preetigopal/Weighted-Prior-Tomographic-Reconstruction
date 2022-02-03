function X_final = AtRadon(b_data,idx,dim,numAngles)
    

% This code is part of the following work which has been submitted to IEEE Transactions on Computational Imaging for peer review.

% Title: "Eliminating object prior-bias from sparse-projection tomographic reconstructions"
% Authors: Preeti Gopal, Sharat Chandran, Imants Svalbe and Ajit Rajwade 

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
X = dct2(backProjImg);
X = reshape(X,[size(backProjImg,1)*size(backProjImg,2) 1]);       


X_final = X;






