%this script compares sparse coding between applying the filter 
%and ways to avoid it. 

% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
num_dict = size(D,3);

imagename = 'lena';

% Load test image
s = single(rgbtogrey(stdimage(imagename)))/255; %use barbara or lena here
s = imresize(s,.5);

% Highpass filter test image

npd = 16;
fltlmbd = 5;
[sl, sh] = lowpass(s, fltlmbd, npd);
shn = sh;
opt = [];
I = randi(255,2,ceil(256*256*.08));
opt.L2Weights = ones(256,256);
for i = 1:size(I,2)
    opt.L2Weights(I(1,i),I(2,i)) = 0;
    shn(I(1,i),I(2,i)) = 1;
end


opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.RelaxParam_z = 1.6;
opt.RelaxParam_z = 1.2;
opt.AutoSigma = 1;
lambda = 0.01; 



[X,~] = cbpdn_l2weight(D,shn,lambda,opt);


scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

figure;
imagesc(sh);
title('original');

figure;
imagesc(shn);
title('noisy');

figure;
imagesc(scnv(D,X));
title('l2reweighed');

