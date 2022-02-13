function store_video(data,fileName,outDirectory,dataset)

data = data - min(data(:));
data = data./max(data(:));
dataPermuted1 = permute(data,[3 1 2]);
dataPermuted2 = permute(data,[3 2 1]);
name = sprintf('%s/%s_%s.mat',outDirectory,dataset,fileName);
save(name,'data');
name = sprintf('%s/%s_%s_z.avi',outDirectory,dataset,fileName);
writevideo (name, data,8);
name = sprintf('%s/%s_%s_y.avi',outDirectory,dataset,fileName);
writevideo (name, dataPermuted1, 8);
name = sprintf('%s/%s_%s_x.avi',outDirectory,dataset,fileName);
writevideo (name, dataPermuted2, 8);



