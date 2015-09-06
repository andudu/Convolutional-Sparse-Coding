function [Y,el_out] = lasso_fista(L, Z, lambda,mu,rho, opt)
%  Solve:        argmin_{Y} (mu/2)Y'LY + (rho/2) ||Y-Z||^2
%                                   +lambda \sum_k ||Y||_1
% via Fista. Spatial coordinates should be flattened before using.
% L is any matrix, preferably sparse
% Z: Load vector (m_k x n_k) x m
% lambda, mu, rho: penalty parameters
% opt.l2w is the Laplacian weights

if(~isempty(opt.Y))
    Y = opt.Y;
else
    Y = zeros(size(Z));
end

if(~isempty(opt.eta))
    eta = opt.eta;
else
    eta = 1.2;
end

if(isfield(opt,'el'))
    el = opt.el;
else
    el = .5;
end

if(isfield(opt,'LinesearchPeriod'))
    lp = opt.LinesearchPeriod;
else
    lp = 5;
end


if(~isfield(opt,'L1Weight'))
    opt.L1Weight = 1;
end

if(~isfield(opt,'tol'))
    opt.tol = 1e-4;
end

Xprv = Y;
flag = 0;
%displaying
if(opt.verbose == 1)
    hstr = 'Itn   Fnc       DFid      l1        r        el    ';
    sfms = '%4d %9.2e %9.2e %9.2e %9.2e %9.2e';
    nsep = 54;
    disp(hstr);
    disp(char('-' * ones(1,nsep)));
    flag = 1;
end


%main iteration

t = 1; k = 0;
r = inf;


while k < opt.MaxMainIter && r>opt.tol;
    linesrch = 1;
    Y_scaled = bsxfun(@times, Y,opt.l2w); 
    Dyz = Y-Z;
    if (mod(k,lp) == 0)
        Fy = .5*mu*trace(Y'*L*Y_scaled)+.5*rho*sum(Dyz(:).^2);
    end
    foo1 = mu*L*Y_scaled;
    foo2 = rho*Dyz;
    G = foo1+foo2;
    while linesrch,
        V = Y- 1/el*G;  %forward gradient step
        X = shrink(V, (lambda/el)*opt.L1Weight);
        if (mod(k,lp) == 0)
            X_scaled = bsxfun(@times, X,opt.l2w); 
            Dxz = X-Z;
            F = .5*mu*trace(X'*L*X_scaled)+ .5*rho*sum(Dxz(:).^2);
            Dxy = X - Y;
            Q =  Fy + Dxy(:)'*G(:) + (el/2)*sum(Dxy(:).^2);
            if F <= Q | el >= 1e5,
                linesrch = 0;
            else
                el = eta*el;
            end
        else
            linesrch = 0; %pretend its already done
        end
    end
    tprv = t;
    t = (1 + sqrt(1 + 4*tprv^2))/2;
    Y = X + ((tprv - 1)/t)*(X - Xprv);
    k = k+1;
    %displaying
    if flag,
        Jl1 = sum(abs(vec(opt.L1Weight .* X)));
        Jf = F + lambda*Jl1;
        r = norm(X(:) - Xprv(:))/norm(X(:));
        disp(sprintf(sfms, k, Jf, F, Jl1, r,el));
    end
    
    %updating to next iteration
    r = norm(X(:) - Xprv(:))/norm(X(:));
    Xprv = X;
end

el_out = el/eta;

end


function u = shrink(v, lambda)

  if isscalar(lambda),
    u = sign(v).*max(0, abs(v) - lambda);
  else
    u = sign(v).*max(0, bsxfun(@minus, abs(v), lambda));
  end

return
end


