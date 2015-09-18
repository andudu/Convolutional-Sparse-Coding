% find a good set of 6 images to learn a nice dictionary

% load the image
ind = [5,16,25,38,46];
load('FlickrCC_512_512.mat');
S0 = S(:,:,ind); 
S0 = single(S0)/255;
clear S; 


% Filter input images and compute highpass images
npd = 16;
fltlmbd = 4;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);

% Construct initial dictionary
numdict = 25;
D0 = zeros(12,12,numdict, 'single');
D0(4:9,4:9,:) = single(randn(6,6,numdict));


% Set up cbpdndliu parameters
lambda = 0.2;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 100;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(Sh,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.6;
opt.DRelaxParam = 1.6;

% Do dictionary learning
[D, X, optinf] = cbpdndliu(D0, Sh, lambda, opt);


% Display learned dictionary
figure;
imdisp(tiledict(D));


