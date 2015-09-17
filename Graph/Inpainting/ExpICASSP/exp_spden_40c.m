addpath('../SPORCO');
sporco
addpath('../TIP_Experiments/LibAAIR');
aairpath

%s0 = single(stdimage('lena.grey'));
%s0 = s0 / 255;
%sn = imnoise(s0,'salt & pepper', 0.3);
load('spnoise30.mat');

% 16x16x128
C = load('../TIP_Experiments/ExpDictLearn2/exprmnt07b.mat');
di = zeros(16,16); di(1,1) = 1;
ds = zeros(16,16); ds(1,1) = 1;
D = cat(3, di, ds, C.D);
clear C;

lambda = 5.13e-03;
mu = 8.2e-02;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 1000;
opt.RelStopTol = 1e-3;
opt.rho = 50*lambda + 1;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;
opt.RelaxParam = 1.8;
opt.L1Weight = single(ones(1, 1, size(D,3), size(sn,3)));
opt.L1Weight(:,:,1,:) = 3.79e-01;

[X, Z, optinf] = cbpdnlc(D, sn, lambda, mu, opt);
DX = ifft2(fft2(D, size(X,1), size(X,2)) .* fft2(X), 'symmetric');
sd = sum(DX(:,:,2:end), 3) + Z;
essim = ssim(s0, sd);
epsnr = psnr(s0, sd);
esnr = snr(s0, sd);
disp(sprintf('SNR: %.2fdB  PSNR: %.2fdB  SSIM: %.3f', esnr, epsnr, essim));

% SNR: 18.95dB  PSNR: 32.20dB  SSIM: 0.915


oi = rmfield(optinf, {'X', 'Xf', 'Y', 'U'});
save('exp_spden_40c.mat', 's0', 'sn', 'sd', 'essim', 'epsnr', 'esnr', ...
     'lambda', 'mu', 'Z', 'X', 'oi');


exit
