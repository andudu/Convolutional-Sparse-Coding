function [D, Y, Y_bar, optinf] = joint_multilearn(D0, S,S_bar, lambda, opt)

% varsplit_multilearn -- Convolutional BPDN Dictionary Learning (Interleaved Update)
%   2 sets of levels at this point
%         argmin_{x_k,d_k} (1/2)||\sum_k h_k * x_k - s||_2^2 +
%                           lambda \sum_k ||x_k||_1
%
%         The solution is computed using Augmented Lagrangian methods
%         (see boyd-2010-distributed) with efficient solution of the 
%         main linear systems (see wohlberg-2014-efficient).
%
% Usage:
%       [D, X, optinf] = cbpdndliu(D0, S, lambda, opt)
%
% Input:
%       D0          Initial dictionary for coarse&fine level
%       S           Input images of larger size
%       S_bar       Input images of smaller size
%       lambda      Regularization parameter
%       opt         Options/algorithm parameters structure (see below)
%
% Output:
%       D           Multires Dictionary filter set (3D array)
%       Y           Larger Coefficient maps (3D array)
%       Y_bar       Smaller Coefficient maps (3D array)
%       optinf      Details of optimisation
%
%
% Options structure fields:
%   HighFilterNum    Number of dictionaries exclusively for high filters
%   Verbose          Flag determining whether iteration status is
%   displayed.
%                    Fields are iteration number, functional value,
%                    data fidelity term, l1 regularisation term, and
%                    primal and dual residuals (see Sec. 3.3 of
%                    boyd-2010-distributed). The values of rho and sigma
%                    are also displayed if options request that they are
%                    automatically adjusted.
%   MaxMainIter      Maximum main iterations
%   AbsStopTol       Absolute convergence tolerance (see Sec. 3.3.1 of
%                    boyd-2010-distributed)
%   RelStopTol       Relative convergence tolerance (see Sec. 3.3.1 of
%                    boyd-2010-distributed)
%   L1Weight         Weight array for L1 norm
%   ImSize           Specifying the actual size of the image. 
%   Y0               Initial value for Y
%   U0               Initial value for U
%   G0               Initial value for G (overrides D0 if specified)
%   H0               Initial value for H
%   rho              Augmented Lagrangian penalty parameter
%   AutoRho          Flag determining whether rho is automatically updated
%                    (see Sec. 3.4.1 of boyd-2010-distributed)
%   AutoRhoPeriod    Iteration period on which rho is updated
%   RhoRsdlRatio     Primal/dual residual ratio in rho update test
%   RhoScaling       Multiplier applied to rho when updated
%   AutoRhoScaling   Flag determining whether RhoScaling value is
%                    adaptively determined (see wohlberg-2015-adaptive). If
%                    enabled, RhoScaling specifies a maximum allowed
%                    multiplier instead of a fixed multiplier
%   sigma            Augmented Lagrangian penalty parameter
%   AutoSigma        Flag determining whether sigma is automatically
%                    updated (see Sec. 3.4.1 of boyd-2010-distributed)
%   AutoSigmaPeriod  Iteration period on which sigma is updated
%   SigmaRsdlRatio   Primal/dual residual ratio in sigma update test
%   SigmaScaling     Multiplier applied to sigma when updated
%   AutoSigmaScaling Flag determining whether SigmaScaling value is
%                    adaptively determined (see wohlberg-2015-adaptive). If
%                    enabled, SigmaScaling specifies a maximum allowed
%                    multiplier instead of a fixed multiplier.
%   StdResiduals     Flag determining whether standard residual definitions
%                    (see Sec 3.3 of boyd-2010-distributed) are used instead
%                    of normalised residuals (see wohlberg-2015-adaptive)
%   XRelaxParam      Relaxation parameter (see Sec. 3.4.3 of
%                    boyd-2010-distributed) for X update
%   DRelaxParam      Relaxation parameter (see Sec. 3.4.3 of
%                    boyd-2010-distributed) for D update
%   LinSolve         Linear solver for main problem: 'SM' or 'CG'
%   MaxCGIter        Maximum CG iterations when using CG solver
%   CGTol            CG tolerance when using CG solver
%   CGTolAuto        Flag determining use of automatic CG tolerance
%   CGTolFactor      Factor by which primal residual is divided to obtain CG
%                    tolerance, when automatic tolerance is active
%   NoBndryCross     Flag indicating whether all solution coefficients
%                    corresponding to filters crossing the image boundary
%                    should be forced to zero.
%   DictFilterSizes  Array of size 2 x M where each column specifies the
%                    filter size (rows x columns) of the corresponding
%                    dictionary filter
%   ZeroMean         Force learned dictionary entries to be zero-mean
%
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2015-04-09
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'Copyright' and 'License' files
% distributed with the library.


if nargin < 5,
  opt = [];
end
checkopt(opt, defaultopts([]));
opt = defaultopts(opt);




% Set up status display for verbose operation
hstr = ['Itn    Fnc      l1       l1_b      '...
        'rx       sx      rxb       sxb      rd      sd       rdb      sdb       rl       sl      '];
sfms = '%4d %5.2e %5.2e %5.2e %5.2e %5.2e %5.2e %5.2e %5.2e %5.2e %5.2e %5.2e %5.2e %5.2e ';
nsep = 145;
if opt.AutoRho,
  hstr = [hstr ' rho     rhob     '];
  sfms = [sfms ' %5.2e %5.2e '];
  nsep = nsep + 10;
end
if opt.AutoSigma,
  hstr = [hstr 'sig       sigb'];
  sfms = [sfms ' %5.2e %5.2e'];
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
if size(S,3) > 1,
  xsz = [size(S,1) size(S,2) size(D0,3) size(S,3)];
  xsz_bar = [size(S_bar,1) size(S_bar,2) size(D0,3)-opt.HighFilterNum size(S_bar,3)];
  hrm = [1 1 1 size(S,3)];
  % Insert singleton 3rd dimension (for number of filters) so that
  % 4th dimension is number of images in input s volume
  S = reshape(S, [size(S,1) size(S,2) 1 size(S,3)]);
  S_bar = reshape(S_bar, [size(S_bar,1) size(S_bar,2) 1 size(S_bar,3)]);  
else
  xsz = [size(S,1) size(S,2) size(D0,3) 1];
  xsz_bar = [size(S_bar,1) size(S_bar,2) size(D0,3)-opt.HighFilterNum 1];
  hrm = 1;
end
xrm = [1 1 size(D0,3)];
Nx = prod(xsz);
Nx_bar = prod(xsz_bar);
Nd = prod(xsz(1:2))*size(D0,3);
Nd_bar = prod(xsz_bar(1:2))*size(D0,3);
cgt = opt.CGTol;

% Dictionary size may be specified when learning multiscale
% dictionary
if isempty(opt.DictFilterSizes),
  dsz = [size(D0,1) size(D0,2)];
else
  dsz = opt.DictFilterSizes;
end

% Mean removal and normalisation projections
Pzmn = @(x) bsxfun(@minus, x, mean(mean(x,1),2));
Pnrm = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(x.^2, 1), 2)));

% Projection of filter to full image size and its transpose
% (zero-pad and crop respectively)
Pzp = @(x) zpad(x, xsz(1:2));
PzpT = @(x) bcrop(x, dsz);
Pzp_bar = @(x) zpad(x, xsz_bar(1:2));

% Projection of dictionary filters onto constraint set
if opt.ZeroMean,
  Pcn = @(x) Pnrm(Pzp(Pzmn(PzpT(x))));
else
  Pcn = @(x) Pnrm(Pzp(PzpT(x)));
end

if opt.ZeroMean,
  Pcn_bar = @(x) Pnrm(Pzp_bar(Pzmn(PzpT(x))));
else
  Pcn_bar = @(x) Pnrm(Pzp_bar(PzpT(x)));
end


% Start timer
tstart = tic;

% Project initial dictionary onto constraint set
D = Pnrm(D0);

% Compute signal in DFT domain
Sf = fft2(S);
Sf_bar = fft2(S_bar);

% Set up algorithm parameters and initialise variables
rho = opt.rho;
delta = opt.delta;
if isempty(rho), rho = 50*lambda+1; end;
rho_bar = rho;


if opt.AutoRho,
  asgr = opt.RhoRsdlRatio;
  asgm = opt.RhoScaling;
end
sigma = opt.sigma;
if isempty(sigma), sigma = size(S,3); end;
sigma_bar = sigma;
if opt.AutoSigma,
  asdr = opt.SigmaRsdlRatio;
  asdm = opt.SigmaScaling;
end
optinf = struct('itstat', [], 'opt', opt);
rx = Inf;
sx = Inf;
rd = Inf;
sd = Inf;
eprix = 0;
eduax = 0;
eprid = 0;
eduad = 0;
rx_bar = Inf;
sx_bar = Inf;
rd_bar = Inf;
sd_bar = Inf;
eprix_bar = 0;
eduax_bar = 0;
eprid_bar = 0;
eduad_bar = 0;



% Initialise main working variables
X = [];
X_bar = [];

%additional dual variable L
L = zeros(size(D,1),size(D,2),size(D,3)-opt.HighFilterNum);

if isempty(opt.Y0), %initialize Y
  Y = zeros(xsz, class(S));
else
  Y = opt.Y0;
end
Yprv = Y;
if isempty(opt.Y0_bar),
  Y_bar = zeros(xsz_bar, class(S_bar));
else
  Y_bar = opt.Y0_bar;
end
Yprv_bar = Y_bar;


if isempty(opt.U0), %initialize U
  if isempty(opt.Y0),
    U = zeros(xsz, class(S));
  else
    U = (lambda/rho)*sign(Y);
  end
else
  U = opt.U0;
end
if isempty(opt.U0_bar),
  if isempty(opt.Y0_bar),
    U_bar = zeros(xsz_bar, class(S_bar));
  else
    U_bar = (lambda/rho_bar)*sign(Y_bar);
  end
else
  U_bar = opt.U0_bar;
end


Df_bar = []; %G_bar
if isempty(opt.G0_bar),
  G_bar = Pzp_bar(D(:,:,1:end-opt.HighFilterNum));
else
  G_bar = opt.G0_bar;
end
Gprv_bar = G_bar;
if isempty(opt.H0_bar),
  if isempty(opt.G0_bar),
    H_bar = zeros(size(G_bar), class(S_bar));
  else
    H_bar = G_bar;
  end
else
  H_bar = opt.H0_bar;
end
Gf_bar = fft2(G_bar, size(S_bar,1), size(S_bar,2));
GSf_bar = bsxfun(@times, conj(Gf_bar), Sf_bar);

Df = []; %G
if isempty(opt.G0),
  G = Pzp(D);
else
  G = opt.G0;
end
Gprv = G;
if isempty(opt.H0),
  if isempty(opt.G0),
    H = zeros(size(G), class(S));
  else
    H = G;
  end
else
  H = opt.H0;
end
Gf = fft2(G, size(S,1), size(S,2));
GSf = bsxfun(@times, conj(Gf), Sf);




% Main loop
k = 1;
while k <= opt.MaxMainIter & ((rx > eprix|sx > eduax|rd > eprid|sd >eduad)|...
   (rx_bar > eprix_bar|sx_bar > eduax_bar|rd_bar > eprid_bar|sd_bar >eduad_bar))     

  % Solve X subproblem. It would be simpler and more efficient (since the
  % DFT is already available) to solve for X using the main dictionary
  % variable D as the dictionary, but this appears to be unstable. Instead,
  % use the projected dictionary variable G
  Xf = solvedbi_sm(Gf, rho, GSf + rho*fft2(Y - U));
  X = ifft2(Xf, 'symmetric');
  clear Xf Gf GSf;

  Xf_bar = solvedbi_sm(Gf_bar, rho_bar, GSf_bar + rho_bar*fft2(Y_bar - U_bar));
  X_bar = ifft2(Xf_bar, 'symmetric');
  clear Xf_bar Gf_bar GSf_bar;  
  
  % See pg. 21 of boyd-2010-distributed
  if opt.XRelaxParam == 1,
    Xr = X;
    Xr_bar = X_bar;
  else
    Xr = opt.XRelaxParam*X + (1-opt.XRelaxParam)*Y;
    Xr_bar = opt.XRelaxParam*X_bar + (1-opt.XRelaxParam)*Y_bar;
  end

  % Solve Y subproblem
  Y = shrink(Xr + U, (lambda/rho)*opt.L1Weight);
  Y_bar = shrink(Xr_bar + U_bar, (lambda/rho_bar)*opt.L1Weight); %set fixed for now
  %crop the Y's according to size. 
  
  if opt.NoBndryCross, %ignore this stuff for now
    Y((end-size(D,1)+2):end,:,:,:) = 0;
    Y(:,(end-size(D,1)+2):end,:,:) = 0;
  end
  
  Yf = fft2(Y);
  Yf_bar = fft2(Y_bar);
  YSf = sum(bsxfun(@times, conj(Yf), Sf), 4);
  YSf_bar = sum(bsxfun(@times, conj(Yf_bar), Sf_bar), 4);
  
  % Update dual variable corresponding to X, Y
  U = U + Xr - Y;
  U_bar = U_bar+Xr_bar-Y_bar;
  clear Xr Xr_bar;

  % Compute primal and dual residuals and stopping thresholds for X update
  nX = norm(X(:)); nY = norm(Y(:)); nU = norm(U(:));
  nX_bar = norm(X_bar(:)); nY_bar = norm(Y_bar(:)); nU_bar = norm(U_bar(:));
  if opt.StdResiduals,
    % See pp. 19-20 of boyd-2010-distributed
    rx = norm(vec(X - Y));
    sx = norm(vec(rho*(Yprv - Y)));
    eprix = sqrt(Nx)*opt.AbsStopTol+max(nX,nY)*opt.RelStopTol;
    eduax = sqrt(Nx)*opt.AbsStopTol+rho*nU*opt.RelStopTol;
    rx_bar = norm(vec(X_bar - Y_bar));
    sx_bar = norm(vec(rho_bar*(Yprv_bar - Y_bar)));
    eprix_bar = sqrt(Nx_bar)*opt.AbsStopTol+max(nX_bar,nY_bar)*opt.RelStopTol;
    eduax_bar = sqrt(Nx_bar)*opt.AbsStopTol+rho_bar*nU*opt.RelStopTol;
  else
    % See wohlberg-2015-adaptive
    rx = norm(vec(X - Y))/max(nX,nY);
    sx = norm(vec(Yprv - Y))/nU;
    eprix = sqrt(Nx)*opt.AbsStopTol/max(nX,nY)+opt.RelStopTol;
    eduax = sqrt(Nx)*opt.AbsStopTol/(rho*nU)+opt.RelStopTol;
    rx_bar = norm(vec(X_bar - Y_bar))/max(nX_bar,nY_bar);
    sx_bar = norm(vec(Yprv_bar - Y_bar))/nU_bar;
    eprix_bar = sqrt(Nx_bar)*opt.AbsStopTol/max(nX_bar,nY_bar)+opt.RelStopTol;
    eduax_bar = sqrt(Nx_bar)*opt.AbsStopTol/(rho*nU)+opt.RelStopTol;
  end
  clear X;

  % Compute l1 norm of Y
  Jl1 = sum(abs(vec(opt.L1Weight .* Y)));
  Jl1_bar = sum(abs(vec(opt.L1Weight .* Y_bar)));
  
  % Update record of previous step Y
  Yprv = Y;
  Yprv_bar = Y_bar;

  
  % Solve D subproblem. Similarly, it would be simpler and more efficient to
  % solve for D using the main coefficient variable X as the coefficients,
  % but it appears to be more stable to use the shrunk coefficient variable Y
  if strcmp(opt.LinSolve, 'SM'),
    Df = solvemdbi_ism(Yf, sigma, YSf + sigma*fft2(G - H));
    Df_bar = solvemdbi_ism(Yf_bar, sigma_bar, YSf_bar + sigma_bar*fft2(G_bar - H_bar));
  else
    [Df, cgst] = solvemdbi_cg(Yf, sigma_bar, YSf + sigma_bar*fft2(G - H), ...
                              cgt, opt.MaxCGIter, Df(:)); %ignore this too
  end
  %clear YSf YSf_bar;
  D = ifft2(Df, 'symmetric');
  D_bar = ifft2(Df_bar,'symmetric');
  if strcmp(opt.LinSolve, 'SM'), clear Df Df_bar; end

  % See pg. 21 of boyd-2010-distributed
  if opt.DRelaxParam == 1,
    Dr = D;
    Dr_bar = D_bar;
  else
    Dr = opt.DRelaxParam*D + (1-opt.DRelaxParam)*G;
    Dr_bar = opt.DRelaxParam*D_bar + (1-opt.DRelaxParam)*G_bar;
  end

  % Solve G G_bar subproblem
  temp = Dr + H; a = temp(:,:,1:end-opt.HighFilterNum); b = temp(:,:,end-opt.HighFilterNum+1:end);
  G(:,:,1:end-opt.HighFilterNum) = Pzp(Pnrm((PzpT(sigma*(a))+delta*(PzpT(G_bar)-L))/(delta+sigma)));    
  G(:,:,end-opt.HighFilterNum+1:end) = Pcn(b);
  clear temp a b;
  Gf = fft2(G);
  GSf = bsxfun(@times, conj(Gf), Sf);
  
  G_bar = Pzp_bar(Pnrm((PzpT(sigma_bar*(Dr_bar + H_bar))+delta*(PzpT(G(:,:,1:end-opt.HighFilterNum))+L))/(delta+sigma)));  
  Gf_bar = fft2(G_bar);
  GSf_bar = bsxfun(@times, conj(Gf_bar), Sf_bar);
  
  % Update dual variable corresponding to D, G
  H = H + Dr - G;
  H_bar = H_bar + Dr_bar - G_bar;
  temp = PzpT(G);
  L = L+temp(:,:,1:end-opt.HighFilterNum)-PzpT(G_bar);
  clear Dr Dr_bar temp;

  % Compute primal and dual residuals and stopping thresholds for D update
  nD = norm(D(:)); nG = norm(G(:)); nH = norm(H(:));  
  nD_bar = norm(D_bar(:)); nG_bar = norm(G_bar(:)); nH_bar = norm(H_bar(:));
  nL = norm(L(:));
  if opt.StdResiduals,
    % See pp. 19-20 of boyd-2010-distributed
    rd = norm(vec(D - G));
    sd = norm(vec(sigma*(Gprv - G)));
    eprid = sqrt(Nd)*opt.AbsStopTol+max(nD,nG)*opt.RelStopTol;
    eduad = sqrt(Nd)*opt.AbsStopTol+sigma*nH*opt.RelStopTol;
    rd_bar = norm(vec(D_bar - G_bar));
    sd_bar = norm(vec(sigma_bar*(Gprv_bar - G_bar)));
    rl = norm(vec(PzpT(G(:,:,1:end-opt.HighFilterNum))-PzpT(G_bar)));
    sl = norm(vec(Gprv_bar-G_bar));   
    eprid_bar = sqrt(Nd_bar)*opt.AbsStopTol+max(nD_bar,nG_bar)*opt.RelStopTol;
    eduad_bar = sqrt(Nd_bar)*opt.AbsStopTol+sigma_bar*nH_bar*opt.RelStopTol;
  else
    % See wohlberg-2015-adaptive
    rd = norm(vec(D - G))/max(nD,nG);
    sd = norm(vec(Gprv - G))/nH;
    eprid = sqrt(Nd)*opt.AbsStopTol/max(nD,nG)+opt.RelStopTol;
    eduad = sqrt(Nd)*opt.AbsStopTol/(sigma_bar*nH)+opt.RelStopTol;
    rd_bar = norm(vec(D_bar - G_bar))/max(nD_bar,nG_bar);
    sd_bar = norm(vec(Gprv_bar - G_bar))/nH_bar;
    rl = norm(vec(PzpT(G(:,:,1:end-opt.HighFilterNum))-PzpT(G_bar)))/max(nG,nG_bar);
    sl = norm(vec(Gprv_bar-G_bar))/nL;
    eprid_bar = sqrt(Nd_bar)*opt.AbsStopTol/max(nD_bar,nG_bar)+opt.RelStopTol;
    eduad_bar = sqrt(Nd_bar)*opt.AbsStopTol/(sigma_bar*nH_bar)+opt.RelStopTol;
  end

%   %max pooling the residuals(how to do that properly? This is not good)
%   rd = sqrt(rd^2+rd_bar^2);
%   sd = sqrt(sd^2+sd_bar^2);
%   eprid = min(eprid,eprid_bar);
%   eduad = min(eduad,eduad_bar);
%   rx = max(rx,rx_bar);
%   sx = max(sx,sx_bar);
%   eprix = min(eprix,eprix_bar);
%   eduax = min(eduax,eduax_bar);  
  
 
  % Apply CG auto tolerance policy if enabled
  if opt.CGTolAuto && (rd/opt.CGTolFactor) < cgt,
    cgt = rd/opt.CGTolFactor;
  end

  clear D;

  % Update record of previous step G
  Gprv = G;
  Gprv_bar = G_bar;


  % Compute data fidelity term in Fourier domain (note normalisation)
  Jdf = sum(vec(abs(sum(bsxfun(@times,Gf,Yf),3)-Sf).^2))/(2*xsz(1)*xsz(2));
  Jdf_bar = sum(vec(abs(sum(bsxfun(@times,Gf_bar,Yf_bar),3)-Sf_bar).^2))/(2*xsz_bar(1)*xsz_bar(2));
  clear Yf Yf_bar;
  Jfn = Jdf +Jdf_bar+ lambda*(Jl1+Jl1_bar);


  % Record and display iteration details
  tk = toc(tstart);
  optinf.itstat = [optinf.itstat;...
       [k Jfn Jl1 rx sx rd sd eprix eduax eprid eduad rho sigma tk]];
  if opt.Verbose,
    dvc = [k,Jfn,Jl1,Jl1_bar, rx, sx,rx_bar,sx_bar, rd, sd,rd_bar,sd_bar,rl,sl];
    if opt.AutoRho,
      dvc = [dvc rho, rho_bar];
    end
    if opt.AutoSigma,
      dvc = [dvc sigma, sigma_bar];
    end
    disp(sprintf(sfms, dvc));
  end

  % See wohlberg-2015-adaptive and pp. 20-21 of boyd-2010-distributed
  if opt.AutoRho,
    if k ~= 1 && mod(k, opt.AutoRhoPeriod) == 0,
      if opt.AutoRhoScaling,
        rhomlt = sqrt(rx/sx);
        if rhomlt < 1, rhomlt = 1/rhomlt; end
        if rhomlt > opt.RhoScaling, rhomlt = opt.RhoScaling; end
      else
        rhomlt = opt.RhoScaling;
      end
      rsf = 1;
      if rx > opt.RhoRsdlRatio*sx, rsf = rhomlt; end
      if sx > opt.RhoRsdlRatio*rx, rsf = 1/rhomlt; end
      rho = rsf*rho;
      U = U/rsf; 
    end
  end
  
  if opt.AutoRho_bar, % for bars too
    if k ~= 1 && mod(k, opt.AutoRhoPeriod) == 0,
      if opt.AutoRhoScaling_bar,
        rhomlt_bar = sqrt(rx_bar/sx_bar);
        if rhomlt_bar < 1, rhomlt_bar = 1/rhomlt_bar; end
        if rhomlt_bar > opt.RhoScaling_bar, rhomlt_bar = opt.RhoScaling_bar; end
      else
        rhomlt_bar = opt.RhoScaling_bar;
      end
      rsf_bar = 1;
      if rx_bar > opt.RhoRsdlRatio*sx, rsf_bar = rhomlt_bar; end
      if sx_bar > opt.RhoRsdlRatio*rx, rsf_bar = 1/rhomlt_bar; end
      rho_bar = rsf_bar*rho_bar;
      U_bar = U_bar/rsf_bar; %for bars too
    end
  end  
  
  if opt.AutoSigma,
    if k ~= 1 && mod(k, opt.AutoSigmaPeriod) == 0,
      if opt.AutoSigmaScaling,
        sigmlt = sqrt(rd/sd);
        if sigmlt < 1, sigmlt = 1/sigmlt; end
        if sigmlt > opt.SigmaScaling, sigmlt = opt.SigmaScaling; end
      else
        sigmlt = opt.SigmaScaling;
      end
      ssf = 1;
      if rd > opt.SigmaRsdlRatio*sd, ssf = sigmlt; end
      if sd > opt.SigmaRsdlRatio*rd, ssf = 1/sigmlt; end
      sigma = ssf*sigma;
      H = H/ssf;
    end
  end

  if opt.AutoSigma_bar,
    if k ~= 1 && mod(k, opt.AutoSigmaPeriod_bar) == 0,
      if opt.AutoSigmaScaling_bar,
        sigmlt_bar = sqrt(rd_bar/sd_bar);
        if sigmlt_bar < 1, sigmlt_bar = 1/sigmlt_bar; end
        if sigmlt_bar > opt.SigmaScaling_bar, sigmlt_bar = opt.SigmaScaling_bar; end
      else
        sigmlt_bar = opt.SigmaScaling_bar;
      end
      ssf_bar = 1;
      if rd_bar > opt.SigmaRsdlRatio*sd_bar, ssf_bar = sigmlt_bar; end
      if sd_bar > opt.SigmaRsdlRatio*rd_bar, ssf_bar = 1/sigmlt_bar; end
      sigma_bar = ssf_bar*sigma_bar;
      H_bar = H_bar/ssf_bar;
    end
  end

  if opt.AutoDelta,
      if k ~= 1 && mod(k, opt.AutoDeltaPeriod) == 0,
          if opt.AutoDeltaScaling,
              sigmlt_d = sqrt(rl/sl);
              if sigmlt_d < 1, sigmlt_d = 1/sigmlt_d; end
              if sigmlt_d > opt.DeltaScaling, sigmlt_d = opt.DeltaScaling; end
          else
              sigmlt_d = opt.DeltaScaling;
          end
          ssf = 1;
          if rl > opt.DeltaRsdlRatio*sl, ssf = sigmlt_d; end
          if sl > opt.DeltaRsdlRatio*rl, ssf = 1/sigmlt_d; end
          delta = ssf*delta;
          L = L/ssf;
      end
  end  
  
  k = k + 1;

end

D = PzpT(G);

% Record run time and working variables
optinf.runtime = toc(tstart);
optinf.Y = Y;
optinf.U = U;
optinf.G = G;
optinf.H = H;
optinf.lambda = lambda;
optinf.rho = rho;
optinf.sigma = sigma;
optinf.cgt = cgt;
if exist('cgst'), optinf.cgst = cgst; end

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


function u = zpad(v, sz)

  u = zeros(sz(1), sz(2), size(v,3), size(v,4), class(v));
  u(1:size(v,1), 1:size(v,2),:,:) = v;

return


function u = bcrop(v, sz)

  if numel(sz) <= 2,
    if numel(sz) == 1
      cs = [sz sz];
    else
      cs = sz;
    end
    u = v(1:cs(1), 1:cs(2), :);
  else
    if size(sz,1) < size(sz,2), sz = sz'; end
    cs = max(sz);
    u = zeros(cs(1), cs(2), size(v,3), class(v));
    for k = 1:size(v,3),
      u(1:sz(k,1), 1:sz(k,2), k) = v(1:sz(k,1), 1:sz(k,2), k);
    end
  end

return


function opt = defaultopts(opt)

  if ~isfield(opt,'Verbose'),
    opt.Verbose = 0;
  end
  if ~isfield(opt,'HighFilterNum'),
      opt.HighFilterNum = 0;
  end
  if ~isfield(opt,'MaxMainIter'),
    opt.MaxMainIter = 1000;
  end
  if ~isfield(opt,'AbsStopTol'),
    opt.AbsStopTol = 1e-6;
  end
  if ~isfield(opt,'RelStopTol'),
    opt.RelStopTol = 1e-4;
  end
  if ~isfield(opt,'L1Weight'),
    opt.L1Weight = 1;
  end
  if ~isfield(opt,'Y0'),
    opt.Y0 = [];
  end
  if ~isfield(opt,'U0'),
    opt.U0 = [];
  end
  if ~isfield(opt,'G0'),
    opt.G0 = [];
  end
  if ~isfield(opt,'H0'),
    opt.H0 = [];
  end
  if ~isfield(opt,'Y0_bar'),
    opt.Y0_bar = [];
  end
  if ~isfield(opt,'U0_bar'),
    opt.U0_bar = [];
  end
  if ~isfield(opt,'G0_bar'),
    opt.G0_bar = [];
  end
  if ~isfield(opt,'H0_bar'),
    opt.H0_bar = [];
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
  if ~isfield(opt,'AutoRhoScaling'),
    opt.AutoRhoScaling = 0;
  end
  if ~isfield(opt,'rho_bar'),
    opt.rho_bar = [];
  end
  if ~isfield(opt,'AutoRho_bar'),
    opt.AutoRho_bar = 0;
  end
  if ~isfield(opt,'AutoRhoPeriod_bar'),
    opt.AutoRhoPeriod_bar = 10;
  end
  if ~isfield(opt,'RhoRsdlRatio_bar'),
    opt.RhoRsdlRatio_bar = 10;
  end
  if ~isfield(opt,'RhoScaling_bar'),
    opt.RhoScaling_bar = 2;
  end
  if ~isfield(opt,'AutoRhoScaling_bar'),
    opt.AutoRhoScaling_bar = 0;
  end
  
  
  if ~isfield(opt,'sigma'),
    opt.sigma = [];
  end
  if ~isfield(opt,'AutoSigma'),
    opt.AutoSigma = 0;
  end
  if ~isfield(opt,'AutoSigmaPeriod'),
    opt.AutoSigmaPeriod = 10;
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
  
  if ~isfield(opt,'sigma_bar'),
    opt.sigma_bar = [];
  end
  if ~isfield(opt,'AutoSigma_bar'),
    opt.AutoSigma_bar = 0;
  end
  if ~isfield(opt,'AutoSigmaPeriod_bar'),
    opt.AutoSigmaPeriod_bar = 10;
  end
  if ~isfield(opt,'SigmaRsdlRatio_bar'),
    opt.SigmaRsdlRatio_bar = 10;
  end
  if ~isfield(opt,'SigmaScaling_bar'),
    opt.SigmaScaling_bar = 2;
  end
  if ~isfield(opt,'AutoSigmaScaling_bar'),
    opt.AutoSigmaScaling_bar = 0;
  end  
  
  if ~isfield(opt,'delta'),
    opt.delta = .5;
  end
  if ~isfield(opt,'AutoDelta'),
    opt.AutoDelta = 0;
  end
  if ~isfield(opt,'AutoDeltaPeriod'),
    opt.AutoDeltaPeriod = 5;
  end
  if ~isfield(opt,'DeltaRsdlRatio'),
    opt.DeltaRsdlRatio = 10;
  end
  if ~isfield(opt,'DeltaScaling'),
    opt.DeltaScaling = 2;
  end
  if ~isfield(opt,'AutoDeltaScaling'),
    opt.AutoDeltaScaling = 0;
  end   
  
  if ~isfield(opt,'StdResiduals'),
    opt.StdResiduals = 0;
  end
  if ~isfield(opt,'XRelaxParam'),
    opt.XRelaxParam = 1;
  end
  if ~isfield(opt,'DRelaxParam'),
    opt.DRelaxParam = 1;
  end
  if ~isfield(opt,'LinSolve'),
    opt.LinSolve = 'SM';
  end
  if ~isfield(opt,'MaxCGIter'),
    opt.MaxCGIter = 1000;
  end
  if ~isfield(opt,'CGTol'),
    opt.CGTol = 1e-3;
  end
  if ~isfield(opt,'CGTolAuto'),
    opt.CGTolAuto = 0;
  end
  if ~isfield(opt,'CGTolAutoFactor'),
    opt.CGTolFactor = 50;
  end
  if ~isfield(opt,'NoBndryCross'),
    opt.NoBndryCross = 0;
  end
  if ~isfield(opt,'DictFilterSizes'),
    opt.DictFilterSizes = [];
  end
  if ~isfield(opt,'ZeroMean'),
    opt.ZeroMean = 0;
  end

return