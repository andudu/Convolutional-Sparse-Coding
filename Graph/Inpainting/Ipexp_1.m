function [ Xnl,Xcn ] = Ipexp_1( D,Sh,L,lambda,mu,ind,maxiter )
%IPEXP_1 Summary of this function goes here
%   Detailed explanation goes here


opt.Verbose = 0;
opt.MaxMainIter = maxiter;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.AutoRhoPeriod = 4;


delta = zeros(size(D,1),size(D,2));
delta(1,1) = 1;
D(:,:,end+1) = delta;
numdict = size(D,3);
opt.L1Weight = ones(size(Sh,1),size(Sh,2),numdict);
opt.L1Weight(:,:,numdict) = 5*ones(size(Sh,1),size(Sh,2));
opt.LapWeight = [ones(1,numdict-1),0]; 

for k = 1:size(ind,2)
    opt.L1Weight(ind(1,k),ind(2,k),numdict) = 0;
end

[Xnl,~] = cbpdnL_lasso(D,Sh,L,lambda,mu,opt);

opt = rmfield(opt, 'LapWeight');

[Xcn,~] = cbpdn(D,Sh,lambda,opt); 

end

