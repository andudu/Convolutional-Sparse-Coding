% script for doing a parameter search on inpainting

lambda_all = [0.01,0.005,.0025, 0.0001];
mu_all = [.1, .05,.025]; 
perc_noise = [.6,.7,.75];
numdict = [15,25,36];

maxiter = 400;

clear D;
% % % Dict 1
% load('CacheData/Dict_12x12x35.mat');
% D = double(D);
% load ([sporco_path,'/Graph/CacheData/Lena_WindowKnearest_varn/Mat/CosineM4.mat']);



% % Dict 3 The real dictionary
% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D = double(D);



S = imresize(double(stdimage('lena.grey'))/255,.5);
[Sl, Sh_ref] = lowpass(S,5,15);


%generate missing pixels


scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
    fft2(x)),3), 'symmetric');


psnr_rec_nl = []; 
snr_rec_nl = []; 


for i = 1:length(perc_noise)
    Sh = Sh_ref;
    ind = [];
    ind(1,:) = randi(size(Sh,1),1,ceil(size(Sh,1)*size(Sh,2)*perc_noise));
    ind(2,:) = randi(size(Sh,2),1,ceil(size(Sh,1)*size(Sh,2)*perc_noise));
    for t = 1:size(ind,2)
        Sh(ind(1,t),ind(2,t)) = 0;
    end
    for j = 1:length(mu_all)
        mu = mu_all(j);
        for k = 1:length(lambda_all)
            lambda = lambda_all(k);
            for l = 1:length(numdict)
                
                disp([num2str(i), num2str(j), num2str(k), num2str(l)]);
                nd = numdict(l);
                D0 = D(:,:,1:nd);
                
                [Xnl,Xcn] = Ipexp_1( D0,Sh,L,lambda, mu,ind,maxiter);
                %reconstruction
                Sh_rec_nl = scnv(D(:,:,1:1:nd),Xnl(:,:,1:1:nd));
                S_rec_nl = Sh_rec_nl + Sl;
                snr_rec_nl(i,j,k,l)= snr(S,S_rec_nl);
                psnr_rec_nl(i,j,k,l) = psnr(S_rec_nl,S);
                
                Sh_rec_cn = scnv(D(:,:,1:1:nd),Xcn(:,:,1:1:nd));
                S_rec_cn = Sh_rec_cn + Sl;
                snr_rec_cn(i,j,k,l)= snr(S,S_rec_cn);
                psnr_rec_cn(i,j,k,l) = psnr(S_rec_cn,S);
            end
        end
    end
end

save('LenaIp3.mat','snr_rec_nl','psnr_rec_nl','snr_rec_cn','psnr_rec_cn','lambda_all','mu_all','perc_noise','numdict','D');


