function [X, optinf] = cbpdntvx(D, S, mu, opt)

% cbpdntvx -- Convolutional Basis Pursuit DeNoising with TV Regularization
%            of each coefficient map instead of l1 regularization
%
%         argmin_{x_k} (1/2)||\sum_k d_k * x_k - s||_2^2 +
%              mu sum_k || sqrt((G_r x_k).^2 + (G_c x_k).^2) ||_2^2 
%
%         The solution is computed using an ADMM approach (see
%         boyd-2010-distributed) with efficient solution of the main
%         linear systems (see wohlberg-2014-efficient).
%
% Usage:
%       [X, optinf] = cbpdntvx(D, S, mu, opt);
%
% Input:
%       D           Dictionary filter set (3D array)
%       S           Input image
%       mu          Regularization parameter (TV)
%       opt         Algorithm parameters structure  
%
% Output:
%       X           Dictionary coefficient map set (3D array)
%       optinf      Details of optimisation
%
%
% Options structure fields:
%   Verbose          Flag determining whether iteration status is displayed.
%                    Fields are iteration number, functional value,
%                    data fidelity term, l1 regularisation term, and
%                    primal and dual residuals (see Sec. 3.3 of
%                    boyd-2010-distributed). The value of rho is also
%                    displayed if options request that it is automatically
%                    adjusted.
%   MaxMainIter      Maximum main iterations
%   AbsStopTol       Absolute convergence tolerance (see Sec. 3.3.1 of
%                    boyd-2010-distributed)
%   RelStopTol       Relative convergence tolerance (see Sec. 3.3.1 of
%                    boyd-2010-distributed)
%   L2Weight         Weighting array for l2 norm of X (default is zero, i.e. 
%                    no l2 norm). Array should have dimensions corresponding
%                    to the non-spatial dimensions of X since spatial weighting
%                    is not possible (i.e. weighting varies only with filter 
%                    and sample index).
%   TVWeight         Weighting array for filters in TV norm of X
%   Zr0              Initial value for Zr
%   Zc0              Initial value for Zc
%   Vr0              Initial value for Vr
%   Vc0              Initial value for Vc
%   sigma            Augmented Lagrangian parameter
%   AutoSigma        Flag determining whether sigma is automatically
%                    updated (see Sec. 3.4.1 of boyd-2010-distributed)
%   AutoSigmaPeriod  Iteration period on which sigma is updated
%   SigmaRsdlRatio   Primal/dual residual ratio in sigma update test
%   SigmaScaling     Multiplier applied to sigma when updated
%   AutoSigmaScaling Flag determining whether SigmaScaling value is
%                    adaptively determined (see wohlberg-2015-adaptive). If
%                    enabled, SigmaScaling specifies a maximum allowed
%                    multiplier instead of a fixed multiplier.
%   SigmaRsdlTarget  Residual ratio targeted by auto sigma update policy.
%   StdResiduals     Flag determining whether standard residual definitions 
%                    (see Sec 3.3 of boyd-2010-distributed) are used instead
%                    of normalised residuals (see wohlberg-2015-adaptive)
%   NonNegCoef       Flag indicating whether solution should be forced to
%                    be non-negative
%   NoBndryCross     Flag indicating whether all solution coefficients
%                    corresponding to filters crossing the image boundary
%                    should be forced to zero.
%   HighMemSolve     Use more memory for a slightly faster solution
%
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2015-08-06
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'Copyright' and 'License' files
% distributed with the library.

if nargin < 4,
  opt = [];
end
checkopt(opt, defaultopts([]));
opt = defaultopts(opt);

% Set up status display for verbose operation
hstr = ['Itn   Fnc       DFid      TV        R(TV)     S(TV)  '];
sfms = '%4d %9.2e %9.2e %9.2e %9.2e %9.2e';
nsep = 54;
if opt.AutoSigma,
  hstr = [hstr '   sigma   '];
  sfms = [sfms ' %9.2e'];
  nsep = nsep + 10;
end
if opt.Verbose && opt.MaxMainIter > 0,
  disp(hstr);
  disp(char('-' * ones(1,nsep)));
end

% Collapsing of trailing singleton dimensions greatly complicates
% handling of both SMV and MMV cases. The simplest approach would be
% if S could always be reshaped to 4d, with dimensions consisting of
% image rows, image cols, a single dimensional placeholder for number
% of filters, and number of measurements, but in the single
% measurement case the third dimension is collapsed so that the array
% is only 3d.
if size(S,3) > 1 && size(S,4) == 1,
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

% Start timer
tstart = tic;

% Compute filters in DFT domain
Df = fft2(D, size(S,1), size(S,2));
% Convolve-sum and its Hermitian transpose
Dop = @(x) sum(bsxfun(@times, Df, x), 3);
DHop = @(x) bsxfun(@times, conj(Df), x);
% Compute signal in DFT domain
Sf = fft2(S);
% S convolved with all filters in DFT domain
DSf = DHop(Sf);

% Set up l2 weight array
if isscalar(opt.L2Weight),
  wl2 = opt.L2Weight;
else
  wl2 = reshape(opt.L2Weight, [1 1 size(opt.L2Weight,1) size(opt.L2Weight,2)]);
end

% Gradient calculations in DFT domain
grv = [-1 1];
Grf = fft2(grv, size(S,1), size(S,2));
gcv = [-1 1]';
Gcf = fft2(gcv, size(S,1), size(S,2));
Gf = conj(Grf).*Grf + conj(Gcf).*Gcf;
% Gradient operators (in DFT domain) and Hermitian transpose
Grop = @(x) bsxfun(@times, Grf, x);
GrHop = @(x) bsxfun(@times, conj(Grf), x);
Gcop = @(x) bsxfun(@times, Gcf, x);
GcHop = @(x) bsxfun(@times, conj(Gcf), x);

% Set up algorithm parameters and initialise variables
sigma = opt.sigma;
if isempty(sigma), sigma = 10*mu; end;
if isempty(opt.SigmaRsdlTarget),
  if opt.StdResiduals,
    opt.SigmaRsdlTarget = 1;
  else
    opt.SigmaRsdlTarget = 1; % Still to be investigated
  end
end
mwr = sigma*Gf + wl2;
if opt.HighMemSolve,
  C = compute_dbd_sm_C(Df, mwr);
else
  C = [];
end
Nx = prod(xsz);
Ns = Nx;
optinf = struct('itstat', [], 'opt', opt);
r = Inf; s = Inf;
epri = 0; edua = 0;

% Initialise main working variables
X = [];
if isempty(opt.Zr0),
  Zr = zeros(xsz);
else
  Zr = opt.Zr0;
end
Zrprv = Zr;
if isempty(opt.Vr0),
  if isempty(opt.Zr0),
    Vr = zeros(xsz);
  else
    Vr = zeros(xsz); % NB: check this
  end
else
  Vr = opt.Vr0;
end
if isempty(opt.Zc0),
  Zc = zeros(xsz);
else
  Zc = opt.Zc0;
end
Zcprv = Zc;
if isempty(opt.Vc0),
  if isempty(opt.Zc0),
    Vc = zeros(xsz);
  else
    Vc = zeros(xsz); % NB: check this
  end
else
  Vc = opt.Vc0;
end


% Main loop
k = 1;
while k <= opt.MaxMainIter && (r > epri | s > edua),

  % Solve X subproblem
  b = DSf + sigma*(GrHop(fft2(Zr - Vr))) + sigma*(GcHop(fft2(Zc - Vc)));
  Xf = solvedbd_sm(Df, mwr, b, C);
  X = ifft2(Xf, 'symmetric');
  xrrs = rrs(DHop(Dop(Xf)) + bsxfun(@times, mwr, Xf), b);

  % Solve Zr, Zc subproblem
  GrX = ifft2(Grop(Xf), 'symmetric');
  GcX = ifft2(Gcop(Xf), 'symmetric');
  [Zr, Zc] = shrinktv(GrX + Vr, GcX + Vc, (mu/sigma)*opt.TVWeight);

  % Update dual variables
  Vr = Vr + GrX - Zr;
  Vc = Vc + GcX - Zc;

  % Compute data fidelity term in Fourier domain (note normalisation)
  Jdf = sum(vec(abs(sum(bsxfun(@times,Df,Xf),3)-Sf).^2))/(2*xsz(1)*xsz(2));
  if isscalar(wl2) && wl2 == 0,
    Jl2X = 0;
  else
    Jl2X = sum(vec(wl2.*sum(sum(X.^2, 1), 2)))/2.0;
  end
  Jtv = sum(vec(bsxfun(@times, sqrt(GrX.^2 + GcX.^2), opt.TVWeight)));
  Jfn = Jdf + Jl2X + mu*Jtv;

  nX = norm(X(:)); nGrX = norm(GrX(:)); nGcX = norm(GcX(:));
  nGrHVr = norm(vec(ifft2(GrHop(fft2(Vr)), 'symmetric')));
  nGcHVc = norm(vec(ifft2(GcHop(fft2(Vc)), 'symmetric')));
  nZr = norm(Zr(:)); nZc = norm(Zc(:)); 
  if opt.StdResiduals,
    % See pp. 19-20 of boyd-2010-distributed
    rr = norm(vec(GrX - Zr));
    sr = norm(vec(sigma*ifft2(GrHop(fft2(Zrprv - Zr)), 'symmetric')));
    rc = norm(vec(GcX - Zc));
    sc = norm(vec(sigma*ifft2(GcHop(fft2(Zcprv - Zc)), 'symmetric')));
    r = 0.5*(rr + rc); s = 0.5*(sr + sc);
    eprir = sqrt(Ns)*opt.AbsStopTol+max(nGrX,nZr)*opt.RelStopTol;
    eduar = sqrt(Nx)*opt.AbsStopTol+sigma*nGrHVr*opt.RelStopTol;
    epric = sqrt(Ns)*opt.AbsStopTol+max(nGcX,nZc)*opt.RelStopTol;
    eduac = sqrt(Nx)*opt.AbsStopTol+sigma*nGcHVc*opt.RelStopTol;
    epri = 0.5*(eprir + eduar); edua = 0.5*(eduar + eduar);
  else
    % See wohlberg-2015-adaptive
    rr = norm(vec(GrX - Zr))/max(nGrX, nZr);
    sr = norm(vec(ifft2(GrHop(fft2(Zrprv - Zr)), 'symmetric')))/nGrHVr;
    rc = norm(vec(GcX - Zc))/max(nGcX, nZc);
    sc = norm(vec(ifft2(GcHop(fft2(Zcprv - Zc)), 'symmetric')))/nGcHVc;
    r = 0.5*(rr + rc); s = 0.5*(sr + sc);
    eprir = sqrt(Ns)*opt.AbsStopTol/ max(nGrX,nZr)+opt.RelStopTol;
    eduar = sqrt(Nx)*opt.AbsStopTol/(sigma*nGrHVr)+opt.RelStopTol;
    epric = sqrt(Ns)*opt.AbsStopTol/ max(nGcX,nZc)+opt.RelStopTol;
    eduac = sqrt(Nx)*opt.AbsStopTol/(sigma*nGcHVc)+opt.RelStopTol;
    epri = 0.5*(eprir + eduar); edua = 0.5*(eduar + eduar);
  end

  
  % Record and display iteration details 
  optinf.itstat = [optinf.itstat; [k Jfn Jdf Jtv r s epri edua sigma xrrs]];
  if opt.Verbose,
    dvc = [k Jfn Jdf Jtv r s];
    if opt.AutoSigma,
      dvc = [dvc sigma];
    end
    disp(sprintf(sfms, dvc));
  end

  % See wohlberg-2015-adaptive and pp. 20-21 of boyd-2010-distributed
  ssf = 1;
  if opt.AutoSigma,
    if k ~= 1 && mod(k, opt.AutoSigmaPeriod) == 0,
      if opt.AutoSigmaScaling,
        sgmmlt = sqrt(r/(s*opt.SigmaRsdlTarget));
        if sgmmlt < 1, sgmmlt = 1/sgmmlt; end
        if sgmmlt > opt.SigmaScaling, sgmmlt = opt.SigmaScaling; end
      else
        sgmmlt = opt.SigmaScaling;
      end
      ssf = 1;
      if r > opt.SigmaRsdlTarget*opt.SigmaRsdlRatio*s, ssf = sgmmlt; end
      if s > (opt.SigmaRsdlRatio/opt.SigmaRsdlTarget)*r, ssf = 1/sgmmlt; end
      sigma = ssf*sigma;
      Vr = Vr/ssf;
      Vc = Vc/ssf;
    end
  end
  if opt.AutoSigma && ssf ~= 1,
    mwr = sigma*Gf + wl2;
    if opt.HighMemSolve,
      C = compute_dbd_sm_C(Df, mwr);
    end
  end
  
  Zrprv = Zr;
  Zcprv = Zc;
  k = k + 1;

end

% Record run time and working variables
optinf.runtime = toc(tstart);
optinf.Xf = Xf;
optinf.Zr = Zr;
optinf.Zc = Zc;
optinf.Vr = Vr;
optinf.Vc = Vc;
optinf.mu = mu;
optinf.sigma = sigma;

if opt.Verbose && opt.MaxMainIter > 0,
  disp(char('-' * ones(1,nsep)));
end

return


function u = vec(v)

  u = v(:);
  
return


function C = compute_dbd_sm_C(Df, mwr)
 
  cn = bsxfun(@rdivide, Df, mwr);
  cd = sum(Df.*bsxfun(@rdivide, conj(Df), mwr), 3) + 1.0;
  C = bsxfun(@rdivide, cn, cd);
  clear cn cd;

return



function u = shrink(v, lambda)
    
  if isscalar(lambda),
    u = sign(v).*max(0, abs(v) - lambda);
  else
    u = sign(v).*max(0, bsxfun(@minus, abs(v), lambda));
  end

return


function [u, v] = shrinktv(x, y, lambda)

  a = sqrt(x.^2 + y.^2);
  if isscalar(lambda),
    b = max(0, a - lambda);
  else
    b = max(0, bsxfun(@minus, a, lambda));
  end
    
  b(a == 0) = 0;
  a(a == 0) = 1;
  b = b./a;
  
  u = b.*x;
  v = b.*y;

return



function opt = defaultopts(opt)

  if ~isfield(opt,'Verbose'),
    opt.Verbose = 0;
  end
  if ~isfield(opt,'MaxMainIter'),
    opt.MaxMainIter = 1000;
  end
  if ~isfield(opt,'AbsStopTol'),
    opt.AbsStopTol = 0;
  end
  if ~isfield(opt,'RelStopTol'),
    opt.RelStopTol = 1e-4;
  end
  if ~isfield(opt,'L2Weight'),
    opt.L2Weight = 0;
  end
  if ~isfield(opt,'TVWeight'),
    opt.TVWeight = 1;
  end
  if ~isfield(opt,'Zr0'),
    opt.Zr0 = [];
  end
  if ~isfield(opt,'Vr0'),
    opt.Vr0 = [];
  end
  if ~isfield(opt,'Zc0'),
    opt.Zc0 = [];
  end
  if ~isfield(opt,'Vc0'),
    opt.Vc0 = [];
  end
  if ~isfield(opt,'sigma'),
    opt.sigma = [];
  end
  if ~isfield(opt,'AutoSigma'),
    opt.AutoSigma = 1;
  end
  if ~isfield(opt,'AutoSigmaPeriod'),
    opt.AutoSigmaPeriod = 1;
  end
  if ~isfield(opt,'SigmaRsdlRatio'),
    opt.SigmaRsdlRatio = 10;
  end
  if ~isfield(opt,'SigmaScaling'),
    opt.SigmaScaling = 2;
  end
  if ~isfield(opt,'AutoSigmaScaling'),
    opt.AutoSigmaScaling = 0;
  end
  if ~isfield(opt,'SigmaRsdlTarget'),
    opt.SigmaRsdlTarget = [];
  end
  if ~isfield(opt,'StdResiduals'),
    opt.StdResiduals = 0;
  end
  if ~isfield(opt,'HighMemSolve'),
    opt.HighMemSolve = 0;
  end

return
