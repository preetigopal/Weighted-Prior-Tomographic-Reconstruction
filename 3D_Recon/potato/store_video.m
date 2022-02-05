function store_video(data,fileName,outDirectory,dataset)

data = data - min(data(:));
data = data./max(data(:));

name = sprintf('%s/%s_%s.mat',outDirectory,dataset,fileName);
save(name,'data');
name = sprintf('%s/%s_%s_z.avi',outDirectory,dataset,fileName);
writevideo (name, data,8);




