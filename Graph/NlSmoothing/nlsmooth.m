function [ sl ] = nlsmooth(s,mu,L,opt)

% Performs a nonlocal smoothing on the signal s, i.e., solving the problem
% below:
% %      |x-s|^2 + \mu<x,Lx> 

% Input:
%   L: Laplacian in the usual format. L{i,j} is the jth window across the
%   ith image. L.M has sparse matrix format, L.phi, L.E for the eigenvector
%   format. 
%   opt:
%       opt.L1: adds an additional L1 Weighting, 1 for on, 0 for off
%       opt.L1Weight: weights for L1


if isfield(L{1,1},'phi')
    lf = 'e';
else
    lf = 'm';
end

if ~isfield(opt,'L1');
    opt.L1 = 0;
end


if ~opt.L1  %just do the nonlocal smoothing
    sl = [];
    for i = 1:size(L,1)
        for j = 1:size(L,2)
            Ltemp = L{i,j};
            I1 = L{i,j}.ind1(1):L{i,j}.ind2(1);
            I2 = L{i,j}.ind1(2):L{i,j}.ind2(2);
            si = reshape(permute(s(I1,I2,i),[2,1]),length(I1)*length(I2),1);
            
            if lf == 'e'
                phi = Ltemp.phi;
                E = Ltemp.E;
                si_c = phi'*si;
                si_par = phi*(si_c);
                si_perp = s - si_par;
                temp_c = bsxfun(@times,1./(mu*E + 1),si_c) ;
                temp_par = phi*temp_c;
                temp_perp = si_perp*1/(mu+1);
                temp = temp_par + temp_perp;
            else
                temp = pcg(speye(size(Ltemp.M))+ mu*Ltemp.M,si);
            end
            temp = reshape(temp,length(I1),length(I2));
            sl(I1,I2,i) = permute(temp,[2,1,3]);
        end
    end
else
    sl = [];
    for i = 1:size(L,1)
        for j = 1:size(L,2)
            Ltemp = L{i,j};
            I1 = L{i,j}.ind1(1):L{i,j}.ind2(1);
            I2 = L{i,j}.ind1(2):L{i,j}.ind2(2);
            si = reshape(permute(s(I1,I2,i),[2,1]),length(I1)*length(I2),1);
            
            if lf == 'e'
                o = {};
                o.verbose = 0;
                o.MaxMainIter = 20;
                temp = eiglasso(Ltemp.phi,Ltemp.E, si,opt.L1Weight,mu,1,o);  
            else
                o = {};
                o.el = .5;
                o.tol = 1e-4;
                o.MaxMainIter = 20;
                temp = lasso_fista(Ltemp.M,si,opt.L1Weight,mu,1,o);
            end
            temp = reshape(temp,length(I1),length(I2));
            sl(I1,I2,i) = permute(temp,[2,1,3]);
        end
    end
        
end





end

