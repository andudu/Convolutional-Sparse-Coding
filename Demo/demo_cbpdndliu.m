% Script demonstrating usage of the cbpdndliu function.
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2015-04-09
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'Copyright' and 'License' files
% distributed with the library.


% % Training images
% S0 = zeros(512, 512, 5, 'single');
% S0(:,:,1) = single(stdimage('lena.grey')) / 255;
% S0(:,:,2) = single(stdimage('barbara.grey')) / 255;
% S0(:,:,3) = single(stdimage('kiel.grey')) / 255;
% S0(:,:,4) = single(rgb2gray(stdimage('mandrill'))) / 255;
% tmp = single(stdimage('man.grey')) / 255;
% S0(:,:,5) = tmp(101:612, 101:612);
% 
% 
% %Reduce images size to speed up demo script
% tmp = zeros(256, 256, 5, 'single');
% for k = 1:size(S0,3),
%   tmp(:,:,k) = imresize(S0(:,:,k), 0.5);
% end
% S0 = tmp;


%just lena
S0= single(stdimage('lena.grey')) / 255;
S0 = imresize(S0,.5);

% Filter input images and compute highpass images
npd = 16;
fltlmbd = 8;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);

% Construct initial dictionary
numdict = 15;
D0 = zeros(12,12,numdict, 'single');
D0(4:9,4:9,:) = single(randn(6,6,numdict));


% Set up cbpdndliu parameters
lambda = 0.1;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 200;
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

% Plot functional value evolution
figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');
