addpath('../SPORCO');
sporco
addpath('../TIP_Experiments/LibAAIR');
aairpath

s0 = single(stdimage('lena.grey'));
s0 = s0 / 255;
sn = s0 + 0.04 * randn(size(s0));

C = load('../TIP_Experiments/ExpDictLearn2/ConvDict.mat');
Dc = C.ConvDict.Dict{17};
Dc = bsxfun(@minus, Dc, mean(mean(Dc,1),2));
g = gauss2d([size(Dc,1) size(Dc,2)], 2); g = g/norm(g(:));
D = cat(3, g, Dc);
clear C

lambda = 2.144e-02;
mu = 3.151e-02;
tvw = 6.888e-02;
l2w = 8.041e-04;

opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 500;
opt.RelStopTol = 2e-3;
opt.rho = 50*lambda + 1;
opt.sigma = 2e1*mu;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;
opt.AutoRhoScaling = 1;
opt.RhoRsdlRatio = 1.2;
opt.RhoScaling = 100;
opt.RhoRsdlTarget = 1;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 1;
opt.AutoSigmaScaling = 1;
opt.SigmaRsdlRatio = 1.2;
opt.SigmaScaling = 100;
opt.SigmaRsdlTarget = 1;
opt.StdResiduals = 0;
opt.L1Weight = ones(1, 1, size(D,3));
opt.L1Weight(1,1,1) = 0;
opt.TVWeight = ones(1, 1, size(D,3));
opt.TVWeight(1,1,1) = tvw/mu;
opt.L2Weight = l2w;

npd = 32;
snp = padarray(sn, [npd npd], 'symmetric', 'both');
  
[X, optinf] = cbpdntv(D, snp, lambda, mu, opt);

DX = ifft2(fft2(D, size(X,1), size(X,2)) .* fft2(X), 'symmetric');
sdp = sum(DX, 3);
sd = sdp((npd+1):(size(sdp,1)-npd), (npd+1):(size(sdp,2)-npd));
essim = ssim(s0, sd);
epsnr = psnr(s0, sd);
esnr = snr(s0, sd);
disp(sprintf('SNR: %.2fdB  PSNR: %.2fdB  SSIM: %.3f', esnr, epsnr, essim));

% SNR: 20.60dB  PSNR: 33.85dB  SSIM: 0.900

oi = rmfield(optinf, {'Xf', 'Y', 'U', 'Zr', 'Zc', 'Vr', 'Vc'});
save('exprmnt58c.mat', 's0', 'sn', 'sd', 'essim', 'epsnr', 'esnr', ...
     'lambda', 'mu', 'tvw', 'l2w', 'X', 'oi');

