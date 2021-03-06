% script for testing 


lambda = .01;
perc_noise = .5;
maxiter = 300;

clear D;

% % % Dict 1
% load('CacheData/Dict_12x12.mat');
% D = D(:,:,1:10);


% % Dict 2
% load ([sporco_path,'/Graph/CacheData/DictLenaCbpdnAll_Window/DictComp326.mat']);
% D = D2;
% %prune the dictionaries
% Aind1 = [];
% for k = 1: size(D,3)
%     if(nnz(X2(:,:,k))>20)
%         Aind1 = [Aind1,k];
%     end
% end
% D = D(:,:,Aind1);



% % Dict 3 The real dictionary
% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D = double(D);



S = double(rgb2gray(stdimage('peppers')))/255;
[Sl, Sh] = lowpass(S,4,15);


%generate missing pixels
ind = [];
ind(1,:) = randi(size(Sh,1),1,ceil(size(Sh,1)*size(Sh,2)*perc_noise));
ind(2,:) = randi(size(Sh,2),1,ceil(size(Sh,1)*size(Sh,2)*perc_noise));

for i = 1:size(ind,2)
    Sh(ind(1,i),ind(2,i)) = 0;
end

opt.Verbose = 1;
opt.MaxMainIter = maxiter;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.AutoRhoPeriod = 4;


scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
    fft2(x)),3), 'symmetric');

delta = zeros(size(D,1),size(D,2));
delta(1,1) = 1;
D(:,:,end+1) = delta;
numdict = size(D,3);
opt.L1Weight = ones(size(S,1),size(S,2),numdict);
opt.L1Weight(:,:,numdict) = 5*ones(size(S,1),size(S,2));


for k = 1:size(ind,2)
    opt.L1Weight(ind(1,k),ind(2,k),numdict) = 0;
end

[X,~] = cbpdn(D,Sh,lambda,opt);
Sh_rec = scnv(D(:,:,1:1:numdict-1),X(:,:,1:1:numdict-1));
S_rec = Sh_rec + Sl;
snr_rec= snr(S,S_rec);
psnr_rec = psnr(S_rec,S);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; 
imagesc(S_rec); 
title(['psnr = ',num2str(psnr_rec), ' srn = ', num2str(snr_rec)]);

figure; 
imagesc(Sh_rec); 
title(['high freq','psnr = ',num2str(psnr_rec), ' srn = ', num2str(snr_rec)]); 


for i = 1:3
    psnrcn = [psnrcn, max(vec(psnr_rec_cn(i,:,:,:)))]; 
    psnrnl = [psnrnl, max(vec(psnr_rec_nl(i,:,:,:)))];
end




