function [X, optinf] = cbpdntv(D, S, lambda, mu, opt)

% cbpdntv -- Convolutional Basis Pursuit DeNoising with TV Regularization
%
%         argmin_{x_k} (1/2)||\sum_k d_k * x_k - s||_2^2 +
%                 lambda \sum_k ||x_k||_1 +
%                 mu || sqrt((G_r \sum_k x_k).^2 + (G_c \sum_k x_k).^2) ||_2^2 
%
%         The solution is computed using the ADMM approach. For
%         details, see S. Boyd, N. Parikh, E. Chu, B. Peleato, and
%         J. Eckstein, "Distributed Optimization and Statistical
%         Learning via the Alternating Direction Method of
%         Multipliers" (doi:10.1561/2200000016).
%
% Usage:
%       [X, optinf] = cbpdntv(D, S, lambda, mu, opt);
%
% Input:
%       D           Dictionary filter set (3D array)
%       S           Input image
%       lambda      Regularization parameter (l1)
%       mu          Regularization parameter (TV)
%       opt         Algorithm parameters structure  
%
% Output:
%       X           Dictionary coefficient map set (3D array)
%       optinf      Details of optimisation
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2014-04-11

if nargin < 5,
  opt = [];
end
opt = defaultopts(opt);

% Collapsing of trailing singleton dimensions greatly complicates
% handling of both SMV and MMV cases. The simplest approach would be
% if S could always be reshaped to 4d, with dimensions consisting of
% image rows, image cols, a single dimensional placeholder for number
% of filters, and number of measurements, but in the single
% measurement case the third dimension is collapsed so that the array
% is only 3d.
if size(S,3) > 1,
  xsz = [size(S,1) size(S,2) size(D,3) size(S,3)];
  hrm = [1 1 1 size(S,3)];
  % Insert singleton 3rd dimension (for number of filters) so that
  % 4th dimension is number of images in input s volume
  S = reshape(S, [size(S,1) size(S,2) 1 size(S,3)]);
else
  xsz = [size(S,1) size(S,2) size(D,3) 1];
  hrm = 1;
end
xrm = [1 1 size(D,3)];
% Compute filters in DFT domain
Df = fft2(D, size(S,1), size(S,2));
grv = [-1 1];
Grf = fft2(grv, size(S,1), size(S,2));
gcv = [-1 1]';
Gcf = fft2(gcv, size(S,1), size(S,2));
if isscalar(opt.TVWeight),
  opt.TVWeight = opt.TVWeight * ones(size(D,3), 1);
end
wgr = reshape(opt.TVWeight, [1 1 length(opt.TVWeight)]);
GrfW = bsxfun(@times, Grf, wgr);
GcfW = bsxfun(@times, Gcf, wgr);
GfW = bsxfun(@times, sqrt(conj(Grf).*Grf + conj(Gcf).*Gcf), wgr);
% Convolve-sum and its Hermitian transpose
Dop = @(x) sum(bsxfun(@times, Df, x), 3);
DHop = @(x) bsxfun(@times, conj(Df), x);
% Compute signal in DFT domain
Sf = fft2(S);
% S convolved with all filters in DFT domain
DSf = DHop(Sf);


% Default lambda is 1/10 times the lambda value beyond which the
% solution is a zero vector
if nargin < 3 | isempty(lambda),
  b = ifft2(DHop(Sf), 'symmetric');
  lambda = 0.1*max(vec(abs(b)));
end
% Set up algorithm parameters and initialise variables
rho = opt.rho;
if isempty(rho), rho = lambda; end;
if opt.AutoRho,
  asgr = opt.RhoRsdlRatio;
  asgm = opt.RhoScaling;
end
sigma = opt.sigma;
if isempty(sigma), sigma = lambda; end;
if opt.AutoSigma,
  asdr = opt.SigmaRsdlRatio;
  asdm = opt.SigmaScaling;
end

X0 = opt.X0;
Nc = size(D,3);
% Initialise main working variables
if isempty(X0),
  X = [];
  Y = zeros(xsz);
  U = Y;
  Zr = zeros(xsz(1:2));
  Zc = zeros(xsz(1:2));
  Vr = Zr;
  Vc = Zc;
else
  X = X0;
  %Y = shrink(X, lambda/rho);
  Y = X;
  b = ifft2(DHop(Sf), 'symmetric');
  c = ifft2(DHop(Dop(fft2(X))), 'symmetric');
  U = (b + rho*Y - c - rho*X)/rho;
  Zr = zeros(xsz(1:2));
  Zc = zeros(xsz(1:2));
  Vr = Zr;
  Vc = Zc;
end

% Set up status display for verbose operation
hstr = ['Itn   Fnc       DFid      l1        TV         R(l1)     S(l1)', ...
       '    R(TV)     S(TV)  '];
sfms = '%4d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e';
nsep = 84;
if opt.AutoRho | opt.RhoUpdate ~= 1,
  hstr = [hstr '   rho  '];
  sfms = [sfms ' %9.2e'];
  nsep = nsep + 10;
end
if opt.AutoSigma | opt.SigmaUpdate ~= 1,
  hstr = [hstr '   sigma   '];
  sfms = [sfms ' %9.2e'];
  nsep = nsep + 10;
end
if opt.Verbose && opt.MaxMainIter > 0,
  disp(hstr);
  disp(char('-' * ones(1,nsep)));
end

tstart = tic;

optinf = struct('itstat', []);
optinf.opt = opt;
k = 1;
r = Inf;
s = Inf;
epri = 0;
edua = 0;
Yprv = Y;
Zrprv = Zr;
Zcprv = Zc;
while k <= opt.MaxMainIter & (r > epri | s > edua),

  % Solve subproblems and update dual variable
  [X, Xf] = convsolve_sm(rho, Y - U, Df, DSf, sigma, ...
                         Zr - Vr, Zc - Vc, GrfW, GcfW, GfW);
  Y = shrink(X + U, (lambda/rho)*opt.L1Weight);
  if opt.NoBndryCross,
    Y((end-size(D,1)+2):end,:,:,:) = 0;
    Y(:,(end-size(D,1)+2):end,:,:) = 0;
  end
  GrWX = ifft2(sum(bsxfun(@times, GrfW, Xf), 3), 'symmetric');
  GcWX = ifft2(sum(bsxfun(@times, GcfW, Xf), 3), 'symmetric');
  [Zr, Zc] = shrinktv(GrWX + Vr, GcWX + Vc, mu/sigma);

  U = U + X - Y;
  Vr = Vr + GrWX - Zr;
  Vc = Vc + GcWX - Zc;
  
  
  % Compute data fidelity term in Fourier domain (note normalisation)
  Jdf = sum(vec(abs(sum(bsxfun(@times,Df,Xf),3)-Sf).^2))/(2*xsz(1)*xsz(2));
  rg = sum(abs(vec(opt.L1Weight .* X)));
  xfs = sum(bsxfun(@times,wgr,Xf),3);
  Jtv = sum(vec(sqrt(GrWX.^2 + GcWX.^2)));
  J1 = Jdf + lambda*rg + mu*Jtv;
  r = norm(vec(X - Y));
  s = norm(vec(rho*(Yprv - Y)));
  rr = norm(vec(GrWX - Zr));
  sr = norm(vec(sigma*(Zrprv - Zr)));
  rc = norm(vec(GcWX - Zc));
  sc = norm(vec(sigma*(Zcprv - Zc)));
  rtv = (rr + rc)/2;
  stv = (sr + sc)/2;
  
  % See Boyd et al., pp. 19-20
  epri = sqrt(Nc)*opt.AbsStopTol + max(norm(X(:)), norm(Y(:)))*opt.RelStopTol;
  edua = sqrt(Nc)*opt.AbsStopTol + norm(rho*U(:))*opt.RelStopTol;
  % NB: TV splitting variables not yet properly handled here
  
  % Record and display iteration details 
  optinf.itstat = [optinf.itstat; [k J1 Jdf rg Jtv r s rtv stv ...
                      epri edua rho sigma]];
  if opt.Verbose,
    dvc = [k J1 Jdf rg, Jtv r s rtv stv];
    if opt.AutoRho | opt.RhoUpdate ~= 1,,
      dvc = [dvc rho];
    end
    if opt.AutoSigma | opt.SigmaUpdate ~= 1,,
      dvc = [dvc sigma];
    end
    disp(sprintf(sfms, dvc));
  end


  if opt.AutoRho,
    if k ~= 1 && mod(k-1, opt.AutoRhoPeriod) == 0,
      % See Boyd et al., pp. 20-21
      vs = 1;
      if r > asgr*s, vs = asgm; end
      if s > asgr*r, vs = 1/asgm; end
      rho = vs*rho;
      U = U/vs;
    end
  elseif opt.RhoUpdate ~= 1 && opt.RhoUpdate*rho < opt.RhoMax,
    rho = opt.RhoUpdate*rho;
    U = U/opt.RhoUpdate;
  end

  if opt.AutoSigma,
    if  k ~= 1 && mod(k-1, opt.AutoSigmaPeriod) == 0,
      % See Boyd et al., pp. 20-21
      vs = 1;
      if rtv > asdr*stv, vs = asdm; end
      if stv > asdr*rtv, vs = 1/asdm; end
      sigma = vs*sigma;
      Vr = Vr/vs;
      Vc = Vc/vs;
    end
  elseif opt.SigmaUpdate ~= 1 && opt.SigmaUpdate*sigma < opt.SigmaMax,
    sigma = opt.SigmaUpdate*sigma;
    Vr = Vr/opt.SigmaUpdate;
    Vc = Vc/opt.SigmaUpdate;
  end

  
  Yprv = Y;
  Zrprv = Zr;
  Zcprv = Zc;
  k = k + 1;

end

% Record run time and working variables
optinf.runtime = toc(tstart);
optinf.Y = Y;
optinf.U = U;
optinf.lambda = lambda;
optinf.rho = rho;

if opt.Verbose && opt.MaxMainIter > 0,
  disp(char('-' * ones(1,nsep)));
end

return


function u = vec(v)

  u = v(:);
  
return


function u = shrink(v, lambda)
    
  u = sign(v).*max(0, abs(v) - lambda);
  
return


function u = shrink2(v, lambda)
  
  vn = norm(v(:));
  if vn == 0,
    u = 0;
  else
    u = v*(max(0, vn - lambda)/vn);
  end
  
return


function [u, v] = shrinktv(x, y, lambda)

  a = sqrt(x.^2 + y.^2);
  b = max(0, a - lambda);
  
  b(a == 0) = 0;
  a(a == 0) = 1;
  
  u = b.*x./a;
  v = b.*y./a;

return


function [X, Xf] = convsolve_sm(rho, Z, Df, DSf, sigma, Zr, Zc, ...
                                GrfW, GcfW, GfW)

% Solve for X_l in (1/2) || \sum_l d_l * x_l - s ||_2^2 +
%                  (\sigma/2) || g_l * \sum_l w_l x_l ||_2^2 +
%                  (\rho/2) \sum_l || x_l - z_l ||_2^2

  Zf = fft2(Z);
  Zrf = fft2(Zr);
  Zcf = fft2(Zc);
  GrZrf = bsxfun(@times, conj(GrfW), Zrf);
  GcZcf = bsxfun(@times, conj(GcfW), Zcf);
  c = DSf + rho*Zf + sigma*GrZrf + sigma*GcZcf;
  
  cGfW = conj(GfW);
  cDf = conj(Df);
  alpha = 1.0 ./ (rho + sigma*sum(bsxfun(@times, GfW, cGfW), 3));
  bc = sum(bsxfun(@times, GfW, c), 3);
  abc = bsxfun(@times, alpha, bc);
  ba = sum(bsxfun(@times, GfW, cDf), 3);
  aba = bsxfun(@times, alpha, ba);
  beta = (c -  bsxfun(@times, sigma*abc, cGfW)) / rho;
  rho = (cDf - bsxfun(@times, sigma*aba, cGfW)) / rho;
  sigma = 1 + sum(bsxfun(@times, Df, rho), 3);
  abeta = sum(bsxfun(@times, Df, beta), 3);
  epsilon = bsxfun(@times, abeta, rho);
  Xf = beta - bsxfun(@rdivide, epsilon, sigma);

  X = ifft2(Xf, 'symmetric');

return


function opt = defaultopts(opt)

  if ~isfield(opt,'Verbose'),
    opt.Verbose = 0;
  end
  if ~isfield(opt,'MaxMainIter'),
    opt.MaxMainIter = 1000;
  end
  if ~isfield(opt,'AbsStopTol'),
    opt.AbsStopTol = 1e-4;
  end
  if ~isfield(opt,'RelStopTol'),
    opt.RelStopTol = 1e-4;
  end
  if ~isfield(opt,'L1Weight'),
    opt.L1Weight = 1;
  end
  if ~isfield(opt,'TVWeight'),
    opt.TVWeight = 1;
  end
  if ~isfield(opt,'X0'),
    opt.X0 = [];
  end
  if ~isfield(opt,'NoBndryCross'),
    opt.NoBndryCross = 0;
  end
  if ~isfield(opt,'rho'),
    opt.rho = [];
  end
  if ~isfield(opt,'AutoRho'),
    opt.AutoRho = 0;
  end
  if ~isfield(opt,'AutoRhoPeriod'),
    opt.AutoRhoPeriod = 10;
  end
  if ~isfield(opt,'RhoRsdlRatio'),
    opt.RhoRsdlRatio = 10;
  end
  if ~isfield(opt,'RhoScaling'),
    opt.RhoScaling = 2;
  end
  if ~isfield(opt,'RhoMax'),
    opt.RhoMax = inf;
  end
  if ~isfield(opt,'RhoUpdate'),
    opt.RhoUpdate = 1;
  end
  if ~isfield(opt,'sigma'),
    opt.sigma = [];
  end
  if ~isfield(opt,'AutoSigma'),
    opt.AutoSigma = 0;
  end
  if ~isfield(opt,'AutoSigmaPeriod'),
    opt.AutoSigmaPeriod = 5;
  end
  if ~isfield(opt,'SigmaRsdlRatio'),
    opt.SigmaRsdlRatio = 10;
  end
  if ~isfield(opt,'SigmaScaling'),
    opt.SigmaScaling = 2;
  end
  if ~isfield(opt,'SigmaMax'),
    opt.SigmaMax = inf;
  end
  if ~isfield(opt,'SigmaUpdate'),
    opt.SigmaUpdate = 1;
  end

return
