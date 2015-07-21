%temporary script for learning a 12 x 12 dictionary from the flicker database

%Training on flicker images
flicker_ind = [4,5,15,17];
load('Flicker1_512_split.mat');
S0 = [];
k = 0;
for i = flicker_ind
    S0(:,:,(k)*4+1) = S(:,:,(i-1)*4+1);
    S0(:,:,(k)*4+2) = S(:,:,(i-1)*4+2);    
    S0(:,:,(k)*4+3) = S(:,:,(i-1)*4+3);
    S0(:,:,(k)*4+4) = S(:,:,(i-1)*4+4);
    k = k+1;
end
S0 = single(S0)/255;
image_num = length(flicker_ind);
clear S;

%Downsample image
S1 = zeros(128, 128, 4*image_num, 'single');
for k = 1:size(S0,3),
    S1(:,:,k) = imresize(S0(:,:,k),.5);
end

% Filter input images and compute highpass images
npd = 16;
fltlmbd = 5;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);
[Sl1, Sh1] = lowpass(S1, fltlmbd, npd);

% % Construct initial dictionary
numdict = 20;
D0 = zeros(12,12,numdict, 'single');
D0(4:9,4:9,:) = single(randn(6,6,numdict));
D0_bar = zeros(6,6,numdict, 'single');
D0_bar(2:5,2:5,:) = single(randn(4,4,numdict));


% Set up cbpdndliu parameters
lambda = 0.2;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 40;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(Sh,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.5;
opt.DRelaxParam = 1.5;

% Do dictionary learning
[D, X, ~] = cbpdndliu(D0, Sh, lambda, opt);

o1.grey =1;
o1.unifscale =0;
square_plot(D,o1);