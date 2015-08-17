function exprmnt58a()

  addpath('../SPORCO');
  sporco

  s0 = single(stdimage('lena.grey'));
  s0 = s0 / 255;
  sn = s0 + 0.04 * randn(size(s0));

  C = load('../TIP_Experiments/ExpDictLearn2/ConvDict.mat');
  Dc = C.ConvDict.Dict{17};
  Dc = bsxfun(@minus, Dc, mean(mean(Dc,1),2));
  g = gauss2d([size(Dc,1) size(Dc,2)], 2); g = g/norm(g(:));
  D = cat(3, g, Dc);
  clear C

  p0 = [1e-1, 5e-3, 5e-3, 1e-4];
  snrfn = @(x) -gwndensnr(s0, sn, D, x);
  options = optimset('fminsearch');
  options = optimset(options, 'Display', 'off');
  [p,ne] = fminsearch(snrfn, p0, options);
  
  disp(sprintf('lambda: %.3e  mu: %.3e  tvw(1): %.3e  l2w: %.3e  SNR: %.3e',...
               p, -ne));

  save('exprmnt58a.mat', 'p0', 'p', 'ne');

return


function sd = gwnden(D, sn, p)

  lambda = p(1);
  mu = p(2);
  tvw = p(3);
  l2w = p(4);
  
  opt = [];
  opt.Verbose = 0;
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
  
return


function e = gwndensnr(s0, sn, D, p)

  sd = gwnden(D, sn, p);
  e = snr(s0, sd);
  disp(sprintf('%.3e %.3e %.3e %.3e  %.3e', p, e));

return
