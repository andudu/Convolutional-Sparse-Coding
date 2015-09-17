addpath('../SPORCO');
sporco
addpath('../TIP_Experiments/LibAAIR');
aairpath
  
s0 = single(stdimage('lena.grey'));
s0 = s0 / 255;
%s0 = imresize(s0, 0.125);

mg = randn(size(s0));
ms = abs(mg) > 0.5;
%sum(ms(:) == 1)/numel(ms)

sn = s0;
sn(ms == 0) = 0;


% 16x16x192
C = load('../TIP_Experiments/ExpDictLearn2/exprmnt07c.mat');
di = zeros(16,16); di(1,1) = 1;
D = cat(3, di, C.D);
clear C;


lambda = 2e-2; mu = 0.75;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 500;
opt.RelStopTol = 2e-3;
opt.rho = 50*lambda + 1;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;
opt.RhoRsdlTarget = 0.1;
opt.RelaxParam = 1.8;

opt.L1Weight = single(ones(size(sn,1), size(sn,2), size(D,3)));
opt.L1Weight(:,:,1) = 1e1*ms;

[X, Z, optinf] = cbpdnlc(D, sn, lambda, mu, opt);

DX = ifft2(fft2(D, size(X,1), size(X,2)) .* fft2(X), 'symmetric');
s1 = sum(DX(:,:,2:end), 3) + Z;

[snr(s0, s1) snr(s0, sn) snr(s1, sn)]
[ssim(s0, s1) ssim(s0, sn) ssim(sn, s1)]

% 23.0884   -4.7101   -4.8098
%  0.9469    0.0795    0.0752

% sum(vec(X(:,:,1) .* ms))  % should be zero

oi = rmfield(optinf, {'Xf', 'X', 'Y', 'U'});
save('exprmnt60.mat','s0','sn','opt','lambda','mu','D','X','Z','oi');



if 0,

  figure;
  subplot(1,3,1);
  imdisp(s0);
  subplot(1,3,2);
  imdisp(sn);
  subplot(1,3,3);
  imdisp(s1);

  figure;
  imdisp(Z);
  
  figure;
  imdisp(sum(DX(:,:,2:end), 3));
  
  figure;
  imagesc(abs(X(:,:,1)));
  
  figure;
  imagesc(sum(abs(X(:,:,2:end)),3));
  
  figure;
  imagesc(sum(abs(DX(:,:,2:end)),3));
  
  
  
  figure;
  subplot(1,4,1);
  plot(oi.itstat(:,2));
  xlabel('Iterations');
  ylabel('Functional value');
  subplot(1,4,2);
  semilogy(oi.itstat(:,7));
  xlabel('Iterations');
  ylabel('Primal residual');
  subplot(1,4,3);
  semilogy(oi.itstat(:,8));
  xlabel('Iterations');
  ylabel('Dual residual');
  subplot(1,4,4);
  plot(oi.itstat(:,12), 'r');
  hold on;
  plot(oi.itstat(:,14), 'b');
  hold off;
  xlabel('Iterations');
  ylabel('GS sub-iterations');
 
  
end