% Script demonstrating usage of the elnet function.
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2015-03-05
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'Copyright' and 'License' files
% distributed with the library.


% Signal and dictionary size
N = 512;
M = 4*N;
% Number of non-zero coefficients in generator
L = 32;
% Noise level
sigma = 0.5;

% Construct random dictionary and random sparse coefficients
D = randn(N, M);
x0 = zeros(M, 1);
si = randperm(M);
si = si(1:L);
x0(si) = randn(L, 1);
% Construct reference and noisy signal
s0 = D*x0;
s = s0 + sigma*randn(N,1);

% Elastic Net for recovery of sparse representation
lambda = 20;
mu = 2;
opt = [];
opt.Verbose = 1;
opt.rho = 100;
opt.AutoRhoScaling = 0;
opt.RhoScaling = 2;
opt.RhoRsdlRatio = 10;
opt.RelStopTol = 1e-6;
[x1, optinf] = elnet(D, s, lambda, mu, opt);

figure;
plot(x0,'r');
hold on;
plot(x1,'b');
hold off;
legend('Reference', 'Recovered', 'Location', 'SouthEast');
title('Dictionary Coefficients');
