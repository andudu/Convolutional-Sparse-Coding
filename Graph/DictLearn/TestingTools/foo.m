
load('CacheData/Dict_12x12.mat');
D1 = D;
D2 = D + .01*randn(size(D));
Din = {D1,D2};
lambda = .02; 
S = imresize(single(stdimage('lena.grey'))/255,.5);
perc_noise = .5;
maxiter = 100;

[ snr_rec, psnr_rec ] = inptest( Din, lambda, S, perc_noise, maxiter);