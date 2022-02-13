%% Parameter setting %%

% % % % % % Confirm your parameters % % % % % % %
function ParamSetting_sprouts(angStepSize)

param.nx = 130;%number of voxels
param.ny = 130;
param.nz = 130;

param.sx = 25;% mm (real size)
param.sy = 25;% mm 114
param.sz = 24;% mm  

%The real detector panel pixel density (number of pixels)
param.nu = 512;		% number of pixels
param.nv = 512;

% Detector setting (real size)
param.su = 512*0.2*4;	% mm (real size)
param.sv = 512*0.2*4;   % mm

% X-ray source and detector setting
param.DSD = 974.182;    %  Distance source to detector 
param.DSO = 126.644;	%  X-ray source to object axis distance

% angle setting
param.dir = 1;   % gantry rotating direction (clock wise/ counter clockwise)
param.dang = angStepSize; % angular step size (deg)
param.deg =  0:param.dang:359.9; 
param.deg = param.deg*param.dir;
param.nProj = length(param.deg);

param.parker = 1; % data with 360 deg -> param.parker = 0 , data less than 360 deg -> param.parker=1 

% % % % % % Confirm your parameters % % % % % % %
 
% filter='ram-lak','cosine', 'hamming', 'hann' 
param.filter='ram-lak'; % high pass filter

param.dx = param.sx/param.nx; % single voxel size
param.dy = param.sy/param.ny;
param.dz = param.sz/param.nz;
param.du = param.su/param.nu;
param.dv = param.sv/param.nv;

param.off_u = 0; param.off_v = 0; % detector rotation shift (real size)

% % % Geometry calculation % % %
param.xs = [-(param.nx-1)/2:1:(param.nx-1)/2]*param.dx;
param.ys = [-(param.ny-1)/2:1:(param.ny-1)/2]*param.dy;
param.zs = [-(param.nz-1)/2:1:(param.nz-1)/2]*param.dz;

param.us = (-(param.nu-1)/2:1:(param.nu-1)/2)*param.du + param.off_u;
param.vs = (-(param.nv-1)/2:1:(param.nv-1)/2)*param.dv + param.off_v;

param.interptype = 'linear'; % 'linear', 'nearest'

% % % % % % Confirm your parameters % % % % % % %
% Only for Matlab version above 2013b with parallel computing toolbox: Speed dependent on your GPU
% You don't need to install CUDA, but install latest graphics driver.
% only Nvidia GPU cards can use this. otherwise please "param.gpu=0"
% This option is semi-GPU code using built-in Matlab GPU functions: several times faster than CPU
param.gpu = 1;


save('parameters_sprouts.mat','param');















