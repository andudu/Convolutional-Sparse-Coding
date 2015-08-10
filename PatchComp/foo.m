% Testing image inpainting

% Training images 512 x 512 
%s = single(stdimage('lena.grey'))/255 ;

% Lightning image
load('lightning.mat');
s = single(x)/255;
s = s+.1*randn(size(s));

% Highpass filter test image
npd = 16;
fltlmbd = 5;
[sl, sh] = lowpass(s, fltlmbd, npd);


% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('8x8x32');

lambda = .2;
opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 100;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.AutoRhoPeriod = 4;


[X,~] = cbpdn(D,sh,lambda,opt);
sh_rec = convsum(D,X,1:1:32);

imwrite(sh*3,'light_whitenoise.png');
imwrite(sh_rec*4.4,'light_whitenoise_rec.png');
