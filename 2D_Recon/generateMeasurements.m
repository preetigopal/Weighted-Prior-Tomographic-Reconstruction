function [y,idx1,dim] = generateMeasurements(testIm,numAngles,outDirectory)

% This code is part of the following work which has been submitted to IEEE Transactions on Computational Imaging for peer review.

% Title: "Eliminating object prior-bias from sparse-projection tomographic reconstructions"
% Authors: Preeti Gopal, Sharat Chandran, Imants Svalbe and Ajit Rajwade 


% This function generates Radon measurements of the 2D test object for
% a fixed set of views.

fname = sprintf('%s/testIm.mat',outDirectory);
save(fname,'testIm');
name = sprintf('%s/testIm.png',outDirectory);
temp = testIm;
temp = temp - min(temp(:));
temp = temp/max(temp(:));
imwrite(temp,name);
% ------------------------------

dim = size(testIm);

% load('../../../../dir_vectors_3668.mat');
load('dir_vectors_3668.mat');
%%---------------selecting same set of angles everytime
id = 1:numAngles;        
idx = id(1:numAngles);
CompleteAngleSet  = mtt(idx,:);
idx1 = (atan(CompleteAngleSet(:,2)./CompleteAngleSet(:,1)))*180/pi;
%%---------------selecting random set of angles everytime
% idx1 = randi(200,[1 numAngles]);

radProj = radon(testIm,idx1); 
y = reshape(radProj,[size(radProj,1)*size(radProj,2) 1]);
y = y(:);

