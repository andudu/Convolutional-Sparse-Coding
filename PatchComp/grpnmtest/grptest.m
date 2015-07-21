% script for testing group norm

% Load data
s_ref = single(rgbtogrey(stdimage('lena')))/255;
%s_ref = imresize(s_ref,.5);

% Construct initial dictionary
% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D1d = reshape(D,size(D,1)*size(D,2),size(D,3));


[sl,sh] = lowpass(s_ref,5,15);



opt = [];
opt.Verbose = 1;
opt.rho = 100;
opt.MaxMainIter =100;
opt.RelStopTol = 1e-6;
opt.MaxMainIter =150;
opt.RelStopTol = 1e-3;

lambda = .1;
mu = .1;

g = zeros(size(sh));
g(:,1:256) = 1;
g = repmat(g,[1,1,36]);


opt.L1Weight = zeros(size(sh));
opt.L21Weight = zeros(size(sh));
opt.L1Weight(:,257:end) = 1;
opt.L21Weight(:,1:256) = 1;



% Computing 
[X,~] = cbpdngrp(D,sh,lambda,mu,g,opt);



% Visualization
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');
shrec = scnv(D,X);
srec = shrec+sl;

imagesc(srec);

