function [rel_mseVal,rmseval] = rmseVol(resultVol,testVol)

% meanSqError = mean((resultVol - testVol).^2,3);
% meanSqError = mean2(meanSqError);
% msi = mean(testVol.^2,3);
% msi = mean2(msi);
% 
% rmseVal = sqrt(meanSqError)/sqrt(msi);

meanSqError = mean((resultVol(:) - testVol(:)).^2);
msi = mean(testVol(:).^2);

rel_mseVal = sqrt(meanSqError)/sqrt(msi);
rmseval = sqrt(meanSqError);