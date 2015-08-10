% Testing Salt and Pepper denoising
% 
% % Training images 512 x 512 
% s = single(stdimage('lena.grey'))/255 ;

% Lightning image
load('lightning.mat');
s = single(x)/255;

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
opt.L1Weight = ones(size(x,1),size(x,2),37);
opt.L1Weight(:,:,37) = .9*ones(size(x));
[X,~] = cbpdn(D,sh,lambda,opt);
sh_rec = convsum(D,X,1:1:numdict-1);

imwrite(sn,'saltpepperlightning_noisy.png');
imwrite(sh_rec+sl,'saltpepperlightning_rec.png');








