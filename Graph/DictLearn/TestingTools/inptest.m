function [ snr_rec, psnr_rec ] = inptest( Din, lambda, S, perc_noise, maxiter )
% function for calculating the psnr and snr for a fake inpainting problem
% for an array of dictionaries

% Input: 
%   Din: cell array of dictionaries
%   lambda: l1 regularizer for inpainting 
%   S: image; 
%   perc_noise: noise percentage. 

% Output: 
% snr, psnr_rec: recorded snr,psnr


[Sl, Sh] = lowpass(S,4,15);


%generate missing pixels
ind = [];
ind(1,:) = randi(size(Sh,1),1,ceil(size(Sh,1)*size(Sh,2)*perc_noise));
ind(2,:) = randi(size(Sh,2),1,ceil(size(Sh,1)*size(Sh,2)*perc_noise));

for i = 1:size(ind,2)
    Sh(ind(1,i),ind(2,i)) = 0;
end

opt.Verbose = 0;
opt.MaxMainIter = maxiter;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.AutoRhoPeriod = 4;
opt.L1Weight = 1;
snr_rec = [];
psnr_rec = [];


scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');


for i = 1:numel(Din)
    D = Din{i};
    delta = zeros(size(D,1),size(D,2));
    delta(1,1) = 1;
    D(:,:,end+1) = delta;
    numdict = size(D,3);
    opt = rmfield(opt,'L1Weight');
    opt.L1Weight = ones(size(S,1),size(S,2),numdict);
    opt.L1Weight(:,:,numdict) = 5*ones(size(S,1),size(S,2));
    for k = 1:size(ind,2)
        opt.L1Weight(ind(1,k),ind(2,k),numdict) = 0;
    end
    
    [X,~] = cbpdn(D,Sh,lambda,opt);
    Sh_rec = scnv(D(:,:,1:1:numdict-1),X(:,:,1:1:numdict-1));
    S_rec = Sh_rec + Sl;
    snr_rec(i) = snr(S,S_rec);
    psnr_rec(i) = psnr(S_rec,S);
end




end

