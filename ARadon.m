
function b = ARadon(X_data,idx,dim,numAngles)

% This code is part of the following work which has been submitted to IEEE Transactions on Computational Imaging for peer review.

% Title: "Eliminating object prior-bias from sparse-projection tomographic reconstructions"
% Authors: Preeti Gopal, Sharat Chandran, Imants Svalbe and Ajit Rajwade 

startPt = 1;
endPt = dim(1)*dim(2);
X = X_data(startPt:endPt);
X = reshape(X,[dim(1) dim(2)]);
X = idct2(X);

startPt = 1;
endPt = numAngles;
angles = idx(startPt:endPt);

radProj = radon(X,angles);  

b1 = reshape(radProj,[size(radProj,1)*size(radProj,2) 1]);

b = b1;

