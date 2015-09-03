function [D, Y, optinf] = cbpdnLdliu(D0, S, L, mu,lambda, opt)

% cbpdndliu -- Convolutional BPDN Dictionary Learning (Interleaved Update)
%
%         argmin_{x_k,d_k} (1/2)||\sum_k h_k * x_k - s||_2^2 +
%                           lambda \sum_k ||x_k||_1
%
%         The solution is computed using Augmented Lagrangian methods
%         (see boyd-2010-distributed) with efficient solution of the 
%         main linear systems (see wohlberg-2014-efficient).
%
% Usage:
%       [D, X, optinf] = cbpdndLliu(D0, S, L,lambda,mu, opt)
%
% Input:
%       D0          Initial dictionary
%       S           Input images
%       lambda      Regularization parameter
%       mu          Regularization parameter for Laplacian
%       L           The Graph Laplacian Cell structure
%       opt         Options/algorithm parameters structure (see below)
%
% Output:
%       D           Dictionary filter set (3D array)
%       X           Coefficient maps (3D array)
%       optinf      Details of optimisation
%
%
% Options structure fields:
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


if nargin < 6,
  opt = [];
end
checkopt(opt, defaultopts([]));
opt = defaultopts(opt);


if ~isempty(opt.Lformat)
    if strcmp(opt.Lformat, 'Eig')
       lf = 'e';
    else
        lf = 'm';
    end
else
    if isfield(L{1},'phi')
        lf = 'e';
    else
        lf = 'm';
    end    
end

% Set up status display for verbose operation
hstr = ['Itn   Fnc       DFid       l1       JLp     '...
        '  r(X)      s(X)      r(D)      s(D) '];
sfms = '%4d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e';
nsep = 84;

if lf == 'e'
   hstr = [hstr '     r(z)      s(z)  '];
   sfms = [sfms ' %9.2e %9.2e'];
   nsep = nsep + 10;
end

if opt.AutoRho,
  hstr = [hstr '     rho  '];
  sfms = [sfms ' %9.2e'];
  nsep = nsep + 10;
end
if lf == 'e' && opt.AutoRhoBar
  hstr = [hstr '     rhob  '];
  nsep = nsep + 10;
end
if opt.AutoSigma,
  hstr = [hstr '     sigma  '];
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
if size(S,3) > 1,
  xsz = [size(S,1) size(S,2) size(D0,3) size(S,3)];
  hrm = [1 1 1 size(S,3)];
  % Insert singleton 3rd dimension (for number of filters) so that
  % 4th dimension is number of images in input s volume
  S = reshape(S, [size(S,1) size(S,2) 1 size(S,3)]);
else
  xsz = [size(S,1) size(S,2) size(D0,3) 1];
  hrm = 1;
end
xrm = [1 1 size(D0,3)];
Nx = prod(xsz);
Nd = prod(xsz(1:2))*size(D0,3);
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

% Projection of dictionary filters onto constraint set
if opt.ZeroMean,
  Pcn = @(x) Pnrm(Pzp(Pzmn(PzpT(x))));
else
  Pcn = @(x) Pnrm(Pzp(PzpT(x)));
end

% Start timer
tstart = tic;

% Project initial dictionary onto constraint set
D = Pnrm(D0);

% Compute signal in DFT domain
Sf = fft2(S);

% Set up algorithm parameters and initialise variables
rho = opt.rho;
if isempty(rho), rho = 50*lambda+1; end;
if opt.AutoRho,
  asgr = opt.RhoRsdlRatio;
  asgm = opt.RhoScaling;
end
sigma = opt.sigma;
if isempty(sigma), sigma = size(S,3); end;
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

if lf == 'e'
    rz = Inf;
    sz = Inf;
    epriz = 0;
    eduaz = 0;
end


% Initialise main working variables
X = [];
if isempty(opt.Y0),
  Y = zeros(xsz, class(S));
else
  Y = opt.Y0;
end
Yprv = Y;
if isempty(opt.U0),
  if isempty(opt.Y0),
    U = zeros(xsz, class(S));
  else
    U = (lambda/rho)*sign(Y);
  end
else
  U = opt.U0;
end
Df = [];
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


if lf == 'm'
    o = cell(1,length(L)); %cell array of options for mini lasso
    for i = 1:size(L,1)
        for j = 1:size(L,2)
            o{i,j}.MaxMainIter = 10;
            o{i,j}.verbose = 0;
            o{i,j}.Y = [];
            o{i,j}.eta = 1.2;
            o{i,j}.tol = opt.RelStopTol/10;
        end
    end
end

if lf == 'e'    
    if isempty(opt.Z0),
        Z = zeros(xsz);
    else
        Z = opt.Z0;
    end
    Zprv = Z;
    if isempty(opt.V0),
        if isempty(opt.V0),
            V = zeros(xsz);
        else
            V = (lambda/rho)*sign(Y);
        end
    else
        V = opt.V0;
    end    
    rho_bar = rho;   
end



% Main loop
k = 1; flag = 1;
while k <= opt.MaxMainIter & flag,

  % Solve X subproblem. It would be simpler and more efficient (since the
  % DFT is already available) to solve for X using the main dictionary
  % variable D as the dictionary, but this appears to be unstable. Instead,
  % use the projected dictionary variable G
  
  if lf == 'm',
    flag = (rx > eprix|sx > eduax|rd > eprid|sd >eduad);
  else
    flag = (rx > eprix|sx > eduax|rd > eprid|sd >eduad |rz > epriz | sz > eduaz ) ; 
  end
  
  
  if lf == 'm'  % Given Sparse Matrices and use Lasso. Remember to add the CG option here. 
    Xf = solvedbi_sm(Gf, rho, GSf + rho*fft2(Y - U));
  else %Given Eigenvectors and use Splitting
    Xf = solvedbi_sm(Gf, rho+rho_bar, GSf + (rho+rho_bar)*fft2(rho/(rho+rho_bar)*(Y - U)+ rho_bar/(rho+rho_bar)*(Z-V)));  
  end
  
  X = ifft2(Xf, 'symmetric');
  clear Xf Gf GSf;

  % See pg. 21 of boyd-2010-distributed
  if opt.XRelaxParam == 1,
    Xr = X;
  else
    Xr = opt.XRelaxParam*X + (1-opt.XRelaxParam)*Y;
  end

  % Solve Y subproblem
  if lf == 'e'
      Y = shrink(Xr + U, (lambda/rho)*opt.L1Weight);
      if opt.NoBndryCross,
          Y((end-size(D,1)+2):end,:,:,:) = 0;
          Y(:,(end-size(D,1)+2):end,:,:) = 0;
      end
      Yf = fft2(Y);
      YSf = sum(bsxfun(@times, conj(Yf), Sf), 4);
  end

  % Solve Z(orY) subproblem. Either via EigDecomp of MiniLasso
  if lf == 'e'  %reshapes!! Careful with the direction
      a = Xr+V; JL_par = 0; JL_perp = 0; Z = zeros(size(a));
      for i = 1:size(L,1) %across images
          for j = 1:size(L,2) %across image windows
              E = L{i,j}.E;
              phi = L{i,j}.phi;
              I1 = L{i,j}.ind1(1):L{i,j}.ind2(1);
              I2 = L{i,j}.ind1(2):L{i,j}.ind2(2);
              ai = reshape(permute(a(I1,I2,:,i),[2,1,3]),size(a,1)*size(a,2) ...
                  ,size(a,3));
              %solve via eigenvectors
              ai_c = phi'*ai;
              ai_par = phi*(ai_c);
              ai_perp = ai - ai_par;
              temp_c = bsxfun(@times,rho_bar./(mu*E + rho_bar),ai_c) ;
              temp_par = phi*temp_c;
              temp_perp = ai_perp*rho_bar/(mu+rho_bar);
              temp = temp_par + temp_perp;
              temp = reshape(temp,size(a,1),size(a,2),size(a,3));
              Z(I1,I2,:,i) = permute(temp,[2,1,3]);
              if opt.Verbose,
                  JL_par = JL_par + sum(vec(temp_c.*bsxfun(@times,temp_c,E)));
                  JL_perp = JL_perp + sum(vec(temp_perp.^2));
              end
          end
      end
      if opt.Verbose,
          JLp = JL_par + JL_perp;
      end
      clear a temp E phi temp temp_par temp_perp temp_c;
  end
  
  if lf == 'm'  %mini laso
      JLp = 0;
      for i = 1:size(L,1)
          for j = 1:size(L,2)
              Ltemp = L{i,j};
              a = Xr+U;
              I1 = L{i,j}.ind1(1):L{i,j}.ind2(1);
              I2 = L{i,j}.ind1(2):L{i,j}.ind2(2);
              [temp,o{i,j}.el] = lasso_fista(Ltemp.M,reshape(permute(a(I1,I2,:,i),[2,1,3]),size(a,1)*size(a,2) ...
                  ,size(a,3)),lambda,mu,rho,o{i,j}); % warm starting
              if opt.Verbose,
                  JLp = trace(temp'*Ltemp.M * temp) +JLp;
              end
              o{i,j}.Y = temp;
              temp = reshape(temp,size(a,1),size(a,2),size(a,3));
              Y(I1,I2,:,i) = permute(temp,[2,1,3]);
              clear a temp;
          end     
      end
      Yf = fft2(Y);
      YSf = sum(bsxfun(@times, conj(Yf), Sf), 4);
  end
  
  
  % Update dual variable corresponding to X, Y
  U = U + Xr - Y;
  if lf == 'e'
      V = V + Xr - Z;
  end
  
  clear Xr;

  % Compute primal and dual residuals and stopping thresholds for X update
  nX = norm(X(:)); nY = norm(Y(:)); nU = norm(U(:));
  if opt.StdResiduals,
    % See pp. 19-20 of boyd-2010-distributed
    rx = norm(vec(X - Y));
    sx = norm(vec(rho*(Yprv - Y)));
    eprix = sqrt(Nx)*opt.AbsStopTol+max(nX,nY)*opt.RelStopTol;
    eduax = sqrt(Nx)*opt.AbsStopTol+rho*nU*opt.RelStopTol;
  else
    % See wohlberg-2015-adaptive
    rx = norm(vec(X - Y))/max(nX,nY);
    sx = norm(vec(Yprv - Y))/nU;
    eprix = sqrt(Nx)*opt.AbsStopTol/max(nX,nY)+opt.RelStopTol;
    eduax = sqrt(Nx)*opt.AbsStopTol/(rho*nU)+opt.RelStopTol;
  end


  if lf == 'e' %calculate additional residual if eig option
      nZ = norm(Z(:));  nV = norm(V(:));
      if opt.StdResiduals,
          % See pp. 19-20 of boyd-2010-distributed
          rz = norm(vec(X - Z));
          sz = norm(vec(rho_bar*(Zprv - Z)));
          epriz = sqrt(Nx)*opt.AbsStopTol+max(nX,nZ)*opt.RelStopTol;
          eduaz = sqrt(Nx)*opt.AbsStopTol+rho_bar*nV*opt.RelStopTol;
      else
          % See wohlberg-2015-adaptive
          rz = norm(vec(X - Z))/max(nX,nZ);
          sz = norm(vec(Zprv - Z))/nV;
          epriz = sqrt(Nx)*opt.AbsStopTol/max(nX,nZ)+opt.RelStopTol;
          eduaz = sqrt(Nx)*opt.AbsStopTol/(rho_bar*nV)+opt.RelStopTol;
      end
  end

  clear X; 
  
  % Compute l1 norm of Y
  Jl1 = sum(abs(vec(opt.L1Weight .* Y)));

  % Update record of previous step Y
  Yprv = Y;
  if lf == 'e'
      Zprv = Z;
  end


  % Solve D subproblem. Similarly, it would be simpler and more efficient to
  % solve for D using the main coefficient variable X as the coefficients,
  % but it appears to be more stable to use the shrunk coefficient variable Y
  if strcmp(opt.LinSolve, 'SM'),
    Df = solvemdbi_ism(Yf, sigma, YSf + sigma*fft2(G - H)); 
  else
    [Df, cgst] = solvemdbi_cg(Yf, sigma, YSf + sigma*fft2(G - H), ...
                              cgt, opt.MaxCGIter, Df(:));
  end
  clear YSf;
  D = ifft2(Df, 'symmetric');
  if strcmp(opt.LinSolve, 'SM'), clear Df; end

  % See pg. 21 of boyd-2010-distributed
  if opt.DRelaxParam == 1,
    Dr = D;
  else
    Dr = opt.DRelaxParam*D + (1-opt.DRelaxParam)*G;
  end

  % Solve G subproblem
  G = Pcn(Dr + H);
  Gf = fft2(G);
  GSf = bsxfun(@times, conj(Gf), Sf);

  % Update dual variable corresponding to D, G
  H = H + Dr - G;
  clear Dr;

  % Compute primal and dual residuals and stopping thresholds for D update
  nD = norm(D(:)); nG = norm(G(:)); nH = norm(H(:));
  if opt.StdResiduals,
    % See pp. 19-20 of boyd-2010-distributed
    rd = norm(vec(D - G));
    sd = norm(vec(sigma*(Gprv - G)));
    eprid = sqrt(Nd)*opt.AbsStopTol+max(nD,nG)*opt.RelStopTol;
    eduad = sqrt(Nd)*opt.AbsStopTol+sigma*nH*opt.RelStopTol;
  else
    % See wohlberg-2015-adaptive
    rd = norm(vec(D - G))/max(nD,nG);
    sd = norm(vec(Gprv - G))/nH;
    eprid = sqrt(Nd)*opt.AbsStopTol/max(nD,nG)+opt.RelStopTol;
    eduad = sqrt(Nd)*opt.AbsStopTol/(sigma*nH)+opt.RelStopTol;
  end

  % Apply CG auto tolerance policy if enabled
  if opt.CGTolAuto && (rd/opt.CGTolFactor) < cgt,
    cgt = rd/opt.CGTolFactor;
  end
  
  clear D;

  % Update record of previous step G
  Gprv = G;


  % Compute data fidelity term in Fourier domain (note normalisation)
  if opt.Verbose,
      Jdf = sum(vec(abs(sum(bsxfun(@times,Gf,Yf),3)-Sf).^2))/(2*xsz(1)*xsz(2));
      clear Yf;
      Jfn = Jdf + lambda*Jl1 + mu*JLp;
      
      
      % Record and display iteration details
      tk = toc(tstart);
      optinf.itstat = [optinf.itstat;...
          [k Jfn Jdf Jl1 rx sx rd sd eprix eduax eprid eduad rho sigma tk]];
  end
  
   if opt.Verbose,
       dvc = [k, Jfn, Jdf, Jl1, JLp, rx, sx, rd, sd];
       if lf == 'e'
           dvc = [dvc, rz,sz];
       end
       if opt.AutoRho,
           dvc = [dvc rho];
       end
       if lf == 'e'
           if opt.AutoRhoBar
               dvc = [dvc rho_bar];
           end
       end
       if opt.AutoSigma,
           dvc = [dvc sigma];
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

  if(lf == 'e')
      if opt.AutoRhoBar,
          if k ~= 1 && mod(k, opt.AutoRhoBarPeriod) == 0,
              if opt.AutoRhoBarScaling,
                  rhomlt_bar = sqrt(rz/sz);
                  if rhomlt_bar < 1, rhomlt_bar = 1/rhomlt_bar; end
                  if rhomlt_bar > opt.RhoBarScaling, rhomlt_bar = opt.RhoBarScaling; end
              else
                  rhomlt_bar = opt.RhoBarScaling;
              end
              rsf = 1;
              if rz > opt.RhoBarRsdlRatio*sz, rsf = rhomlt_bar; end
              if sz > opt.RhoBarRsdlRatio*rz, rsf = 1/rhomlt_bar; end
              rho_bar = rsf*rho_bar;
              V = V/rsf;
          end
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


% % play movie per iteration
% if k == 1,
%     o1.newfig = 1;
%     o2.newfig = 1;
% else
%     o1.newfig = 0;
%     o2.newfig = 0;
% end
% if k < 50,
%     o1.grey = 1;
%     o1.fighand = square_plot((PzpT(G)),o1);
%     drawnow;
%     o2.sparse = 1;
%     o2.fighand = square_plot((Y),o2);
%     drawnow;
% else
%     if mod(k,5) == 0,
%         o1.grey = 1;
%         o1.fighand = square_plot((PzpT(G)),o1);
%         drawnow;
%         o2.sparse = 1;
%         o2.fighand = square_plot((Y),o2);
%         drawnow;
%     end
% end
k = k + 1;

end

D = PzpT(G);

% Record run time and work(ing variables
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
  if ~isfield(opt,'Z0'),
    opt.Z0 = [];
  end
  if ~isfield(opt,'V0'),
    opt.V0 = [];
  end  
  
  if ~isfield(opt,'G0'),
    opt.G0 = [];
  end
  if ~isfield(opt,'H0'),
    opt.H0 = [];
  end
  if ~isfield(opt,'rho'),
    opt.rho = [];
  end
  if ~isfield(opt,'rho_bar'),
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
  if ~isfield(opt,'Lformat'),
      opt.Lformat = [];
  end
  if ~isfield(opt,'AutoRhoBar'),
    opt.AutoRhoBar = 0;
  end
  if ~isfield(opt,'AutoRhoBarPeriod'),
    opt.AutoRhoBarPeriod = 10;
  end
  if ~isfield(opt,'RhoBarRsdlRatio'),
    opt.RhoBarRsdlRatio = 10;
  end
  if ~isfield(opt,'RhoBarScaling'),
    opt.RhoBarScaling = 2;
  end
  if ~isfield(opt,'AutoRhoBarScaling'),
    opt.AutoRhoBarScaling = 0;
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
