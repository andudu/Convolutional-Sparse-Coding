
% Standard Training images
S0 = double(stdimage('lena.grey'))/255;

% Filter input images and compute highpass images
npd = 16;
fltlmbd = 5;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);




% Construct initial dictionary
load('CacheData/Dict_12x12.mat');
D0 = double(D);
D0 = D0(:,:,1:1:15);
numdict = size(D0,3);



% Cbpdn using a nice dictionary 
lambda = 0.02;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 300;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(Sh,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.XRelaxParam = 1.5;
opt.DRelaxParam = 1.5;
[X1,~, optinf1] = cbpdn(D0,Sh, lambda,opt);


% Set up cbpdndliu parameters
lambda = 0.02;
mu = .2;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;


% Do dictionary learning
[D2, X2, optinf2] = cbpdndliu(D0, Sh, lambda, opt);

a1 = [];
for i = 1:15
    a1(i) = nnz(X1(:,:,i)); 
end
a1 = sort(a1); 

a2 = [];
for i = 1:15
    a2(i) = nnz(X2(:,:,i)); 
end
a2 = sort(a2);   

figure;
plot(1:1:15,a1,'b', 1:1:15,a2,'r'); 





