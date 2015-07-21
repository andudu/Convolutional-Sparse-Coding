% script for testing group norm

% Load data
s_ref = single(rgbtogrey(stdimage('lena')))/255;
s_low = imresize(s_ref,.5);


% Construct initial dictionary
% Load dictionary
load('Flicker20im34dictcold_dict.mat');

[sl_low,sh_low] = lowpass(s_low,10,15);


opt = [];
opt.Verbose = 1;
opt.rho = 150;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 3;
opt.MaxMainIter =100;
opt.RelStopTol = 1e-3;
lambda = .01;


[x2,~] = cbpdn(D_bar,sh_low,lambda,opt);


% Visualization
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

x = PxT(x2); 


                           
shrec = scnv(D,x);
srec_high = 4*shrec+imresize(sl_low,2);

figure;
imagesc(srec_high);
title('superres');
colorbar;

psnr(srec_high,s_ref);



