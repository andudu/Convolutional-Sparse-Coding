function [ Y_bar,V ] = eiglasso( phi, E, Z,lambda,mu, rho,opt)
% solve the Lasso problem via eigenvector decomposition
%         argmin_{Y} (mu/2)Y'LY + (rho/2) ||Y-Z||^2
%                           lambda \sum_k ||Y||_1
% via ADMM. Spatial coordinates should be flattened before using.
% L is symmetric semi-positive definite, having rank r
% phi: (m_k x n_k) x r matrix, top r eigenvector of the matrix L
% E: length r vector,  top r eigenvalues of L
% Z: Load vector (m_k x n_k) x m
% lambda, mu: penalty parameters
% rho: stepsize for ADMM


if(~isempty(opt.Y_bar))
    Y_bar = opt.Y_bar;
else
    Y_bar = zeros(size(Z));
end

if(~isempty(opt.V))
    V = opt.V;
else
    V = zeros(size(Z));
end

flag = 0;
if(opt.verbose == 1)
    hstr = 'Itn   Fnc       DFid      l1        r         s    ';
    sfms = '%4d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e';
    nsep = 54;
    disp(hstr);
    disp(char('-' * ones(1,nsep)));
    flag = 1;
    Y_bar_prv = Y_bar;
end


k = 0;
A = .5*(Z+Y_bar-V);
A_hat = phi'*A;
A_par = phi*A_hat; %project Z onto nontrivial eigenspace
A_perp = A-A_par; %perpendicular part

while (k< opt.Maxiter)
    Y_par_hat = bsxfun(@rdivide,A_hat,(mu*E+2*rho))*2*rho; % Y_par update
    Y = phi*Y_par_hat + (2*rho/(mu+2*rho))*A_perp;
    Y_bar = shrink(Y + V, lambda/rho);
    V = Y-Y_bar+V;
    A = .5*(Z+Y_bar-V);  
    A_hat = phi'*A;
    A_par = phi*A_hat;
    A_perp = A-A_par;
    k = k+1;
    if flag,
        Y_perp = Y-phi*Y_par_hat;
        Jdf = .5*mu* sum(vec(bsxfun(@times,E,Y_par_hat.^2))) + .5*mu* sum(vec(Y_perp.*Y_perp))+.5*rho *sum(vec((Y-Z).^2));
        Jl1 = sum(abs(vec(Y)))*lambda;
        Jfn  = Jl1 + Jdf;
        r = norm(vec(Y - Y_bar));
        s = norm(vec(rho*(Y_bar_prv - Y_bar)));
        disp(sprintf(sfms, k, Jfn, Jdf, Jl1, r, s));    
        Y_bar_prv = Y_bar;
    end
end


end


function u = shrink(v, lambda)

  if isscalar(lambda),
    u = sign(v).*max(0, abs(v) - lambda);
  else
    u = sign(v).*max(0, bsxfun(@minus, abs(v), lambda));
  end

return
end


