
% Standard Training images 512 x 512 lena
S0 = double(stdimage('lena.grey'))/255;

% Filter input images and compute highpass images
npd = 16;
fltlmbd = 5;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);


% % Dict 3 The real dictionary
% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D = double(D);



% Cbpdn using a nice dictionary 
lambda = 0.2;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(Sh,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.XRelaxParam = 1.5;
opt.DRelaxParam = 1.5;
opt.RelStopTol = 1e-4;
[X1,~, optinf1] = cbpdn(D0,Sh, lambda,opt);


% Set up cbpdndliu parameters
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;


% Do dictionary learning
[D2, X2, optinf2] = cbpdndliu(D0, Sh, lambda, opt);

fid_sc = optinf1.itstat(:,3);
fid_dict = optinf2.itstat(:,3); 

l1_sc = optinf1.itstat(:,4);
l1_dict = optinf2.itstat(:,4); 


figure; 
l = min(length(fid_sc),length(fid_dict)); 
plot(1:1:l, fid_sc(1:1:l), 'b',1:1:l, fid_dict(1:1:l), 'r'); 
title(['Fidelity between Sparse coding and dict learning, lambda = ',num2str(lambda)]); 
xlabel('iteration'); 
legend('Sparse Code', 'Dict Learn'); 


figure; 
l = min(length(fid_sc),length(fid_dict)); 
plot(1:1:l, l1_sc(1:1:l), 'b',1:1:l, l1_dict(1:1:l), 'r'); 
title(['L1 between Sparse coding and dict learning, lambda = ',num2str(lambda)]); 
xlabel('iteration'); 
legend('Sparse Code', 'Dict Learn'); 





