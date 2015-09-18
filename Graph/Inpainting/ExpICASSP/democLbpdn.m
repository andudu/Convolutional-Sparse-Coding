% Demo for the new script the effects of a laplacian sparse coding

% script for testing 

% Params
lambda = .01;
mu1 = .75;
mu2 = .1; 
perc_noise = .55;

clear D;



% load Laplacian
load([sporco_path,'/Graph/CacheData/Lena512_Knearest_vark/Mat/CosineM3.mat']);

% % Dict 3 The real dictionary
% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D = double(D);
temp = D;
D(:,:,1) = zeros(12,12); 
D(1,1,1) = 1; 
D(:,:,2:end+1) = temp; 
numdict = size(D,3);

% load image
S = double(stdimage('lena.grey'))/255;
S_ref = S;


%generate missing pixels
ind = [];
ind(1,:) = randi(size(S,1),1,ceil(size(S,1)*size(S,2)*perc_noise));
ind(2,:) = randi(size(S,2),1,ceil(size(S,1)*size(S,2)*perc_noise));

for i = 1:size(ind,2)
    S(ind(1,i),ind(2,i)) = 0;
end






%%%%%%%%%%%%%%%%%%%%%%%% Optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = [];
opt.L1Weight = ones(size(S,1),size(S,2),size(D,3));
opt.L1Weight(:,:,1) = 10*ones(size(S,1),size(S,2));
for t = 1:size(ind,2)
    opt.L1Weight(ind(1,t),ind(2,t),1) = 0;
end
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.RelStopTol = 1e-3;
opt.rho = 50*lambda + 1;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 1;
opt.RhoRsdlTarget = 0.1;
opt.RelaxParam = 1.8;
opt.LapWeight = [0,ones(1,numdict-1)]; 

scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
    fft2(x)),3), 'symmetric');


[Xnl,Znl,~] = cLbpdnlc(D,S,L,lambda,mu1,mu2,opt);
Sh_rec_nl = scnv(D(:,:,2:end),Xnl(:,:,2:end));
S_rec_nl = Sh_rec_nl + Znl;
psnr_rec = psnr(S_rec_nl,S_ref);


[Xcn,Zcn,~] = cbpdnlc(D,S,lambda,mu1,opt); 
Sh_rec_cn = scnv(D(:,:,2:end),Xcn(:,:,2:end));
S_rec_cn = Sh_rec_cn + Zcn;
psnr_rec_cn = psnr(S_rec_cn,S_ref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; 
imagesc(S_rec_nl); 
title(['psnr nl = ',num2str(psnr_rec)]);

figure; 
imagesc(Sh_rec_nl); 
title('high freq nl'); 


figure; 
imagesc(S_rec_cn); 
title(['psnr cn = ',num2str(psnr_rec_cn)]);

figure; 
imagesc(Sh_rec_cn); 
title('high freq cn'); 