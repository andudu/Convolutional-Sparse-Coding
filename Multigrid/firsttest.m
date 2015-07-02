% Script demonstrating usage of the cbpdn function.
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2015-04-09
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'Copyright' and 'License' files
% distributed with the library.


% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
num_dict = 36;
D = dmap('12x12x36');


% Load test image
s = single(rgbtogrey(stdimage('lena')))/255;
if isempty(s),
  error('Data required for demo scripts has not been installed.');
end

% Highpass filter test image
npd = 16;
fltlmbd = 5;



[sl, sh] = lowpass(s, fltlmbd, npd);
sh2 = imresize(sh,0.5);




D2 = zeros(6,6,num_dict);
%dowsize the dictionary
for i = 1:num_dict 
    D2(:,:,2) = imresize(D(:,:,2),0.5);
end

% Compute representation
lambda = 0.03;
opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;

% %multigrid
% tic;
% [X2,U2,optinf]= cbpdn(D2,sh2,lambda,  opt); 
% t2 = toc;
% X0 = zeros(size(X));
% U0 = zeros(size(X));
% for i = 1:num_dict
%     X0(:,:,i) = imresize(X2(:,:,i),2);
%     U0(:,:,i) = imresize(U2(:,:,i),2);
% end
% opt.Y0 = X0/4;
% opt.U0 = U0/4;
% tic;
% opt.MaxMainIter = 20;
% [X,~,optinf] = cbpdn(D,sh,lambda,opt);
% t2 = toc+t2;
% clear X X0  X2 U0 U2 ;

opt.MaxMainIter = 30;


%full stuff
tic;
opt.Y0 = [];
opt.U0 = [];
figure;
image(100*(sh+0.3278));
title('original');
[X,~, optinf] = foocbpdn(D, sh, lambda, opt);


