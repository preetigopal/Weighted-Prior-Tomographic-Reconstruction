function [rel_mseVal,rmseval] = rmseVol(resultVol,testVol)

meanSqError = mean((resultVol(:) - testVol(:)).^2);
msi = mean(testVol(:).^2);

rel_mseVal = sqrt(meanSqError)/sqrt(msi);
rmseval = sqrt(meanSqError);