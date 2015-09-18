% Compare Results for Image Inpainting 256 x 256

% Demo for the new script the effects of a laplacian sparse coding

% script for testing

% Params

noise_level = [.4,.5,.6,.65,.7,.75];
lambda_best = [.005,.01,.01,.01,.02,.02];
mu1_best = [.5,.5,.5,.5,.5,.5];
mu2_all = [.05,.1,.15];

clear D;

psnr_all_nl = [];
psnr_all_cn = [];

for i = 1:length(noise_level)
        disp(num2str(i));

        load(['CorImNoise256',num2str(i),'.mat']);
        mu1 = mu1_best(i);
        % load Laplacian
        load([sporco_path,'/Graph/CacheData/Lena256_Knearest_Ip/Lena256_Knearest_Ip/Mat/CosineM',num2str(i),'.mat']);
        
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
        S_ref = imresize(double(stdimage('lena.grey'))/255,.5);
        S = S_c;
        
        %%%%%%%%%%%%%%%%%%%%%%%% Optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        opt = [];
        opt.L1Weight = ones(size(S,1),size(S,2),size(D,3));
        opt.L1Weight(:,:,1) = 10*ones(size(S,1),size(S,2));
        for t = 1:size(ind,2)
            opt.L1Weight(ind(1,t),ind(2,t),1) = 0;
        end
        opt.Verbose = 0;
        opt.MaxMainIter = 500;
        opt.RelStopTol = 2e-3;
        opt.rho = 50*lambda + 1;
        opt.AutoRho = 1;
        opt.AutoRhoPeriod = 1;
        opt.RhoRsdlTarget = 0.1;
        opt.RelaxParam = 1.8;
        opt.LapWeight = [0,ones(1,numdict-1)];
        
        scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
            fft2(x)),3), 'symmetric');
        
        for j = 1:length(mu2_all)
            mu2 = mu2_all(j);
            [Xnl,Znl,~] = cLbpdnlc(D,S,L,lambda,mu1,mu2,opt);
            Sh_rec_nl = scnv(D(:,:,2:end),Xnl(:,:,2:end));
            S_rec_nl = Sh_rec_nl + Znl;
            psnr_rec_nl = psnr(S_rec_nl,S_ref);
            psnr_all_nl(i,j) = psnr_rec_nl;
        end
        
        
        [Xcn,Zcn,~] = cbpdnlc(D,S,lambda,mu1,opt);
        Sh_rec_cn = scnv(D(:,:,2:end),Xcn(:,:,2:end));
        S_rec_cn = Sh_rec_cn + Zcn;
        psnr_rec_cn = psnr(S_rec_cn,S_ref);
        psnr_all_cn(i) = psnr_rec_cn;
         
end

save('CompareResults256.mat','psnr_all_cn','psnr_all_nl','noise_level','lambda_best','mu1_best','mu2_all');