function fbp = filteredBackProjection2D(y,numAngles,idx1,dim)

% This code is part of the following work which has been submitted to IEEE Transactions on Computational Imaging for peer review.

% Title: "Eliminating object prior-bias from sparse-projection tomographic reconstructions"
% Authors: Preeti Gopal, Sharat Chandran, Imants Svalbe and Ajit Rajwade 


proj = reshape(y,[size(y,1)/numAngles numAngles]);
fbp = iradon(proj,idx1,'linear','Cosine');
fbp = fbp(2:2+dim(1)-1,2:2+dim(2)-1);