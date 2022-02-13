function writeColorVideo (filename, vid, framerate)

[H,W,C,NF] = size(vid);
ff = VideoWriter(filename);
ff.FrameRate = framerate;


open (ff);
for i=1:NF
    writeVideo(ff,vid(:,:,:,i));
end
close (ff);