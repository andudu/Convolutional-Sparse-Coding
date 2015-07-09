% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('8x8x32');
num_dict = size(D,3);

imagename = 'lena';

% Load test image
s = single(rgbtogrey(stdimage(imagename)))/255; %use barbara here
s2 = imresize(s,.5);
s3 = imresize(s,.25);


scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

                                                   
                           
% Highpass filter test image
npd = 16;
fltlmbd = 5;
[sl, sh] = lowpass(s, fltlmbd, npd);
[sl2, sh2] = lowpass(s2, fltlmbd+2, npd);
[sl3, sh3] = lowpass(s3, fltlmbd+4, npd);


%setting up parameters
lambda = 0.01;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 500;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;
[X_d, optinf1] = cbpdn(D, sh3, lambda, opt);

sh_d = scnv(D,X_d);
s_d = sh_d + sl3;
ps_nr_d = psnr(s_d,s3);
figure;
imagesc(s_d); colormap(gray);
title(['psnr of Regular = ',num2str(ps_nr_d)]);

% figure; 
% imagesc(sh);
% figure;
% imagesc(sh2);
% figure;
% imagesc(sh3);
