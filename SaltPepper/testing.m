% Testing Salt and Pepper denoising

% Training images 512 x 512 
s = single(stdimage('lena.grey'))/255 ;


% Add Salt and Pepper
sn = imnoise(s, 'salt & pepper');


% Highpass filter test image
npd = 16;
fltlmbd = 5;
[sl, sh] = lowpass(sn, fltlmbd, npd);


% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
delta = zeros(12,12);
delta(1,1) = 1;
D(:,:,end+1) = delta;
numdict = size(D,3);


lambda = .03;
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.AutoRhoPeriod = 4;


[X,~] = cbpdn(D,sh,lambda,opt);
sh_rec = convsum(D,X,1:1:numdict-1);
% 
% figure;
% imagesc(sn);
% axis off;
% saveas(gcf,'saltpeppernoisylena','png');
% 
% figure;
% imagesc(sh_rec+sl);
% axis off;
% saveas(gcf,'saltpepperreclena','png');

imwrite(sn,'saltpeppernoisylena.png');
imwrite(sh_rec+sl,'saltpepperreclena.png');








