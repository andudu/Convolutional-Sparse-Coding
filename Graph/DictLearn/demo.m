
%generate the eigenvectors using Nystrom

%load([sporco_path,'Graph/CacheData/']);
load([sporco_path,'/Graph/CacheData/Standard/Mat/CosineM1.mat']);
ltemp = L;
L = {};
L{1,1}.M= ltemp{1}.M;
L{1,1}.ind1 = [1,1];
L{1,1}.ind2 = [256,256];

% L{1,1}.phi = phi;
% E = (E./max(E)).^0.75;
% L{1,1}.E = E;
% L{1,1}.ind1 = [1,1];
% L{1,1}.ind2 = [256,256];


% Standard Training images
S0 = [];
S0(:,:,1) = double(imresize(double(stdimage('lena.grey'))/255,.5));

% Filter input images and compute highpass images
npd = 16;
fltlmbd = 5;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);
%Sh = Sh+.06*randn(size(Sh));

% Construct initial dictionary
numdict = 10;
D0 = zeros(12,12,numdict);
D0(4:9,4:9,:) = randn(6,6,numdict);
% load('CacheData/Dict_12x12.mat');
% D0 = double(D);
% D0 = D0(:,:,1:1:numdict);

numdict = size(D0,3);

% Set up cbpdndliu parameters
lambda = 0.2;
mu = .2;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 50;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(Sh,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.5;
opt.DRelaxParam = 1.5;
opt.AutoRhoBar = 1;
opt.AutoRhoBarPeriod = 10;
opt.Lformat = 'Mat';

% Do dictionary learning
[D1, X1, ~] = cbpdnLdliu(D0, Sh,L, lambda,mu, opt);
%[D2, X2, ~] = cbpdndliu(D0, Sh, lambda, opt);

% y = [];
% for i = 1:size(D2,3)
%     y(:,:,i) = convsum(D2,X2,i);
% end
% square_plot(y,{});
% 
% figure; 
% imagesc(convsum(D2,X2,1:1:size(D2,3)));



% % Display learned dictionary
% figure;
% imdisp(tiledict(D1));

% Display learned dictionary

% o.grey = 1;
% square_plot(D1,o);
% p.sparse = 1;
% square_plot(X1,p);

% o.grey = 1;
% square_plot(D2,o);
% p.sparse = 1;
% square_plot(X2,p);

% save(fname,'D1','X1','D2','X2','lambda','mu')