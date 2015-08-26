% generate initial condition for diffusion

%load a saved noise
sref = double(stdimage('lena.grey')) / 255;
sref = imresize(sref,.5);
[sl,sh] = lowpass(sref,5,15);


% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = double(dmap('12x12x36'));


lambda = 0.35;
opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 100;
opt.rho = 10;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;


[X,~] = cbpdn(D,sh,lambda,opt);

X = sum3(abs(X));

figure; imagesc(X);

save('Coef2.mat','X');