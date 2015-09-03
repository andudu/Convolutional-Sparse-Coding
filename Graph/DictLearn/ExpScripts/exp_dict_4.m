% experiment 4 for dictionary learning
% learning from noisy images with nice initialization.
% remember to set D0 to the nice dictionary, 
% maxit = 300, and S0 to be the noisy image for training. 

function [D1,X1, D2, X2 ] = exp_dict_4(D0, S0, L, lambda, mu, maxit)

%%%%%%%%%%%%%%%%%%%%%%%%%% Load the Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load([sporco_path,'Graph/CacheData/']);
% global sporco_path;
% % Standard Training images
% S0 = zeros(512, 512, 5, 'single');
% S0(:,:,1) = single(stdimage('lena.grey')) / 255;
% S0(:,:,2) = single(stdimage('barbara.grey')) / 255;
% S0(:,:,3) = single(stdimage('kiel.grey')) / 255;
% S0(:,:,4) = single(rgb2gray(stdimage('mandrill'))) / 255;
% tmp = single(stdimage('man.grey')) / 255;
% S0(:,:,5) = tmp(101:612, 101:612);

% % Reduce images size to speed up demo script
% tmp = [];
% kk = 1;
% for k = imind,
%   tmp(:,:,kk) = S0(:,:,k);
%   kk = kk+1;
% end
% S0 = tmp;





%%%%%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Filter input images and compute highpass images
npd = 16;
fltlmbd = 4;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);


% Set up cbpdndliu parameters
opt1 = [];
opt1.Verbose = 0;
opt1.MaxMainIter = maxit;
opt1.rho = 50*lambda + 0.5;
opt1.sigma = size(Sh,3);
opt1.AutoRho = 1;
opt1.AutoRhoPeriod = 10;
opt1.AutoSigma = 1;
opt1.AutoSigmaPeriod = 10;
opt1.AutoRhoBar = 1;
opt1.AutoRhoBarPeriod = 10;
opt1.Lformat = 'Sparse';
opt1.XRelaxParam = 1.5;
opt1.DRelaxParam = 1.5;



% Set up Laplacian parameters
opt2 = [];
opt2.Verbose = 0;
opt2.MaxMainIter = maxit;
opt2.rho = 50*lambda + 0.5;
opt2.sigma = size(Sh,3);
opt2.AutoRho = 1;
opt2.AutoRhoPeriod = 10;
opt2.AutoSigma = 1;
opt2.AutoSigmaPeriod = 10;
opt2.XRelaxParam = 1.5;
opt2.DRelaxParam = 1.5;


%%%%%%%%%%%%%%%%%%%%%%%% Dictionary Learning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[D1, X1, ~] = cbpdnLdliu(D0, Sh,L, lambda,mu, opt1);
[D2, X2, ~] = cbpdndliu(D0, Sh, lambda, opt2);

end

