function X_final = AtRadon_wavelet(b_data,idx,dim,numAngles,wname)
    
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
[CA,CH,CV,CD] = dwt2(backProjImg,wname);
CA = imresize(CA,[dim(1)/2 dim(2)/2]);
CH = imresize(CH,[dim(1)/2 dim(2)/2]);
CV = imresize(CV,[dim(1)/2 dim(2)/2]);
CD = imresize(CD,[dim(1)/2 dim(2)/2]);
X = [CA CH;CV CD];

X = reshape(X,[size(backProjImg,1)*size(backProjImg,2) 1]);       


X_final = X;






