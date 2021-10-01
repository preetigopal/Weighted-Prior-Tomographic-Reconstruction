
function b = ARadon_wavelet(X_data,idx,dim,numAngles,wname,size_CA)

startPt = 1;
endPt = dim(1)*dim(2);
X = X_data(startPt:endPt);
X = reshape(X,[dim(1) dim(2)]);
CA_recovered = X(1:size_CA(1),1:size_CA(2));
CH_recovered = X(1:size_CA(1),size_CA(2)+1:end);
CV_recovered = X(size_CA(1)+1:end,1:size_CA(2));
CD_recovered = X(size_CA(1)+1:end,size_CA(2)+1:end);
X = idwt2(CA_recovered,CH_recovered,CV_recovered,CD_recovered,wname);

X = imresize(X,[dim(1) dim(2)]);
startPt = 1;
endPt = numAngles;
angles = idx(startPt:endPt);

radProj = radon(X,angles);  

b1 = reshape(radProj,[size(radProj,1)*size(radProj,2) 1]);

b = b1;

