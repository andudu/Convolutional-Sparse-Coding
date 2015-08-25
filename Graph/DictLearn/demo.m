
%generate the eigenvectors using Nystrom

%load([sporco_path,'Graph/CacheData/']);
load('Euclideaneig1.mat');

S0 = imresize(single(stdimage('lena.grey')) / 255,.5);

% Filter input images and compute highpass images
npd = 16;
fltlmbd = 5;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);
Sh = Sh+.06*randn(size(Sh));

% Construct initial dictionary
numdict = 10;
 D0 = zeros(12,12,numdict, 'single');
 D0(4:9,4:9,:) = single(randn(6,6,numdict));
% load('CacheData/Dict_12x12.mat');
% D0 = double(D);
% D0 = D0(:,:,1:1:numdict);

numdict = size(D0,3);

% Set up cbpdndliu parameters
lambda = 0.2;
mu = .3;
opt = [];
opt.Verbose = 1;
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
E = (E./max(E)).^0.75;
L{1,1}.E = E;
L{1,1}.ind1 = [1,1];
L{1,1}.ind2 = [256,256];
[D1, X1, ~] = cbpdnLdliu(D0, Sh,L, lambda,mu, opt);
[D2, X2, ~] = cbpdndliu(D0, Sh, lambda, opt);

% % Display learned dictionary
% figure;
% imdisp(tiledict(D1));

% Display learned dictionary

o.grey = 1;
square_plot(D1,o);
p.sparse = 1;
square_plot(X1,p);

o.grey = 1;
square_plot(D2,o);
p.sparse = 1;
square_plot(X2,p);

% save(fname,'D1','X1','D2','X2','lambda','mu')