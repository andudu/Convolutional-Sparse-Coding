function [X, optinf] = cbpdnfista(D, S, lambda, opt)

% cbpdn -- Convolutional Basis Pursuit DeNoising
%
%         argmin_{x_k} (1/2)||\sum_k d_k * x_k - s||_2^2 +
%                           lambda \sum_k ||x_k||_1
%
%         The solution is computed using the FISTA approach.
%   
% Usage:
%       [X, optinf] = cbpdnfista(D, S, lambda, opt);
%
% Input:
%       D           Dictionary filter set (3D array)
%       S           Input image
%       lambda      Regularization parameter
%       opt         Algorithm parameters structure  
%
% Output:
%       X           Dictionary coefficient map set (3D array)
%       optinf      Details of optimisation
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2014-05-23

if nargin < 4,
  opt = [];
end
opt = defaultopts(opt);
L = opt.L0;
eta = opt.eta;

% Collapsing of trailing singleton dimensions greatly complicates
% handling of both SMV and MMV cases. The simplest approach would be
% if s could always be reshaped to 4d, with dimensions consisting of
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
  xsz = [size(S,1) size(S,2) size(D,3) size(S,4)];
  hrm = 1;
end
xrm = [1 1 size(D,3)];
% Compute filters in DFT domain
Df = fft2(D, size(S,1), size(S,2));
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

% Initialise main working variables
X0 = opt.X0;
if isempty(X0),
  X = zeros(xsz);
  Xf = zeros(xsz);
  if opt.DFTDomainUpdate,
    Zf = zeros(xsz);
  else
    Z = zeros(xsz);
    Rx = -S;
    Rz = -S;
  end
  Rxf = -Sf;
  Rzf = -Sf;
else
  X = X0;
  Xf = fft2(X);
  Rxf = sum(bsxfun(@times,Df,Xf),3) - Sf;
  Rzf = Rxf;
  if opt.DFTDomainUpdate,
    Zf = Xf;
  else
    Z = X;
    Rx = ifft2(Rxf, 'symmetric');
    Rz = Rx;
  end
end

% Set up status display for verbose operation
hstr = 'Itn   Fnc       DFid      l1        L         t         FrcChng  ';
sfms = '%4d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e';
nsep = 64;
if opt.Verbose && opt.MaxMainIter > 0,
  disp(hstr);
  disp(char('-' * ones(1,nsep)));
end

tstart = tic;

optinf = struct('itstat', []);
optinf.opt = opt;
k = 1;
t = 1;
frcxd = inf;
while k <= opt.MaxMainIter & frcxd > opt.RelStopTol,

  % Gradient
  Gf = bsxfun(@times, conj(Df), Rxf);
  if opt.DFTDomainUpdate,
    Xfprv = Xf;
  else
    G = ifft2(Gf, 'symmetric');
    Xprv = X;
  end
  
  linesrch = 1;
  while linesrch,

    if opt.DFTDomainUpdate,
      Vf = Zf - (1/L)*Gf;
      V = ifft2(Vf, 'symmetric');
      X = shrink(V, (lambda/L)*opt.L1Weight);
      Xf = fft2(X);
      Rxf = sum(bsxfun(@times,Df,Xf),3) - Sf;
      F = sum(abs(Rxf(:)).^2)/2;
      Dxzf = Xf - Zf;
      %Q = sum(abs(Rzf(:)).^2)/2 + ccat(Dxzf(:))'*ccat(Gf(:)) + ...
      %    (L/2)*sum(abs(Dxzf(:)).^2);
      Q = sum(abs(Rzf(:)).^2)/2 + real(Dxzf(:)'*Gf(:)) + ...
          (L/2)*sum(abs(Dxzf(:)).^2);
    else
      V = Z - (1/L)*G;      
      X = shrink(V, (lambda/L)*opt.L1Weight);
      Xf = fft2(X);
      Rxf = sum(bsxfun(@times,Df,Xf),3) - Sf;
      Rx = ifft2(Rxf, 'symmetric');
      F = sum(Rx(:).^2)/2;
      Dxz = X - Z;
      Q = sum(Rz(:).^2)/2 + Dxz(:)'*G(:) + (L/2)*sum(Dxz(:).^2);
    end
    
    if F <= Q,
      linesrch = 0;
    else
      L = eta*L;
    end
      
  end
  
  tprv = t;
  t = (1 + sqrt(1 + 4*tprv^2))/2;
  
  if opt.DFTDomainUpdate,
    Zf = Xf + ((tprv - 1)/t)*(Xf - Xfprv);
  else
    Z = X + ((tprv - 1)/t)*(X - Xprv);
    Zf = fft2(Z);
  end
  Rzf = sum(bsxfun(@times,Df,Zf),3) - Sf;
  
  if ~opt.DFTDomainUpdate,
    Rz = ifft2(Rzf, 'symmetric');
  end
    
  % Compute data fidelity term in Fourier domain (note normalisation)
  Jdf = sum(abs(Rxf(:)).^2)/(2*xsz(1)*xsz(2));
  %Jdf = sum(Rx(:).^2)/2;
  rg = sum(abs(vec(opt.L1Weight .* X)));
  J1 = Jdf + lambda*rg;
  if opt.DFTDomainUpdate,
    frcxd = norm(Xf(:) - Xfprv(:))/norm(Xf(:));
  else
    frcxd = norm(X(:) - Xprv(:))/norm(X(:));
  end
  
  ittm = toc(tstart);
  optinf.itstat = [optinf.itstat; [k J1 Jdf rg L t ittm]];
  if opt.Verbose,
    disp(sprintf(sfms, k, J1, Jdf, rg, L, t, frcxd));
  end

  k = k + 1;
    
end


% Record run time and working variables
optinf.runtime = toc(tstart);
optinf.Gf = Gf;
if opt.DFTDomainUpdate,
  optinf.Zf = Zf;
else
  optinf.Z = Z;
end

if opt.Verbose && opt.MaxMainIter > 0,
  disp(char('-' * ones(1,nsep)));
end



return


function u = vec(v)

  u = v(:);
  
return


function u = ccat(v)

  u = [real(v); imag(v)];
  
return


function u = shrink(v, lambda)
    
  u = sign(v).*max(0, abs(v) - lambda);
  
return

function opt = defaultopts(opt)

  if ~isfield(opt,'Verbose'),
    opt.Verbose = 0;
  end
  if ~isfield(opt,'MaxMainIter'),
    opt.MaxMainIter = 1000;
  end
  if ~isfield(opt,'RelStopTol'),
    opt.RelStopTol = 1e-3;
  end
  if ~isfield(opt,'L1Weight'),
    opt.L1Weight = 1;
  end
  if ~isfield(opt,'X0'),
    opt.X0 = [];
  end
  if ~isfield(opt,'L0'),
    opt.L0 = 1;
  end
  if ~isfield(opt,'eta'),
    opt.eta = 5;
  end
  if ~isfield(opt,'DFTDomainUpdate'),
    opt.DFTDomainUpdate = 1;
  end

return