function [D1,X1,D2,X2] = comptest(D0, lambda,mu)

%generate the eigenvectors using Nystrom

%load([sporco_path,'Graph/CacheData/']);
load('Euclideaneig1.mat');

S0 = imresize(single(stdimage('lena.grey')) / 255,.5);

% Filter input images and compute highpass images
npd = 16;
fltlmbd = 4;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);
Sh = Sh + .05*randn(size(Sh));

% Construct initial dictionary
% numdict = 15;
% D0 = zeros(12,12,numdict, 'single');
% D0(4:9,4:9,:) = single(randn(6,6,numdict));

numdict = size(D0,3);

% Set up cbpdndliu parameters
%lambda = 0.2;
%mu = lambda/2;
opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 300;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(Sh,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.3;
opt.DRelaxParam = 1.3;
opt.AutoRhoBar = 1;
opt.AutoRhoBarPeriod = 10;
opt.Lformat = 'Eig';

% Do dictionary learning
L = {};
L{1,1}.phi = phi;
L{1,1}.E = E;
L{1,1}.ind1 = [1,1];
L{1,1}.ind2 = [size(S0,1),size(S0,2)];
E = (E./max(E)).^0.75;
[D1, X1, ~] = cbpdnLdliu(D0, Sh,L, lambda,mu, opt);
[D2, X2, ~] = cbpdndliu(D0, Sh, lambda, opt);

% % Display learned dictionary
% figure;
% imdisp(tiledict(D1));
% 
% % Display learned dictionary
% figure;
% imdisp(tiledict(D2));

%save(fname,'D1','X1','D2','X2','lambda','mu')