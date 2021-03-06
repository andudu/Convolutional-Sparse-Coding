%experiment 2 for dictionary learning
function [D2, X2,Aind2] = exp_dict_3(D0,lambda, imind,maxit)

%%%%%%%%%%%%%%%%%%%%%%%%%% Load the Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load([sporco_path,'Graph/CacheData/']);
global sporco_path;
% Standard Training images
S0 = zeros(512, 512, 5, 'single');
S0(:,:,1) = single(stdimage('lena.grey')) / 255;
S0(:,:,2) = single(stdimage('barbara.grey')) / 255;
S0(:,:,3) = single(stdimage('kiel.grey')) / 255;
S0(:,:,4) = single(rgb2gray(stdimage('mandrill'))) / 255;
tmp = single(stdimage('man.grey')) / 255;
S0(:,:,5) = tmp(101:612, 101:612);

% Reduce images size to speed up demo script
tmp = [];
kk = 1;
for k = imind,
  tmp(:,:,kk) = S0(:,:,k);
  kk = kk+1;
end
S0 = tmp;





%%%%%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Filter input images and compute highpass images
npd = 16;
fltlmbd = 4;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);
%Sh = Sh + .05*randn(size(Sh));

% Construct initial dictionary
% numdict = 15;
% D0 = zeros(12,12,numdict, 'single');
% D0(4:9,4:9,:) = single(randn(6,6,numdict));
% numdict = size(D0,3);

% Set up cbpdndliu parameters
%lambda = 0.2;
%mu = lambda/2;
opt = [];
opt.Verbose = 0;
opt.MaxMainIter = maxit;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(Sh,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.5;
opt.DRelaxParam = 1.5;


%%%%%%%%%%%%%%%%%%%%%%%% Dictionary Learning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[D2, X2, ~] = cbpdndliu(D0, Sh, lambda, opt);

%%%%%%%%%%%%%%%%%%%%%%%% Recording Unused Dicts %%%%%%%%%%%%%%%%%%%%%%%%%%%
Aind2 = [];
for i = 1: size(D0,3)
    if(nnz(X2(:,:,i,:))>=20)
       Aind2 = [Aind2,i]; 
    end      
end


end

