% Script demonstrating usage of the cbpdndliu function.
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2015-04-09
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'Copyright' and 'License' files
% distributed with the library.

% load the image
ind = [  28 ,   34 ,   42 ,   61 ,   99 ];
load('Flicker1_512_split.mat');
S0 = S(:,:,ind); 
clear S; 


% Filter input images and compute highpass images
npd = 16;
fltlmbd = 8;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);

% Construct initial dictionary
numdict = 25;
D0 = zeros(12,12,numdict, 'single');
D0(4:9,4:9,:) = single(randn(6,6,numdict));


% Set up cbpdndliu parameters
lambda = 0.7;
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

% Plot functional value evolution
figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');
