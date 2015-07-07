%this script generates a noisy image and 

% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
num_dict = size(D,3);

%loda a test image
imagename = 'lena';
s = single(rgbtogrey(stdimage(imagename)))/255; %use barbara or lena here
s = imresize(s, 0.5);
s_ref = s;
%add the noise (20/255 percent peak signal)
sigma = 20/255;
s = s+ sigma*randn(size(s));


%function for convsum
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

% Highpass filter test image
npd = 16;
fltlmbd = 5;
[sl, sh] = lowpass(s, fltlmbd, npd);
[sl_ref,sh_ref] = lowpass(s_ref, fltlmbd, npd);

%setting up params
lambda = 0.15;
mu = 0.3*lambda; %half the size of lambda
lambda2 = 1*(lambda+mu);
opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 400;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;


%parameter plot
% disp('jnt denoising');
% X_nt = cbpdnjnt_new(D,sh,lambda,mu,opt);
% disp('clean image denoising');
% X_dc = cbpdn(D,sh_ref,lambda2,opt);
disp('regular denoising');
X = cbpdn(D,sh,lambda,opt);
% disp('BM3D denoising');
% [ps_nr_bm3d,S_bm3d] = BM3D(s_ref, s);
% disp('BM3D denoising Clean');
% [ps_nr_bm3dc,S_bm3dc] = BM3D(s_ref, s_ref);

%computing 
sh_d = scnv(D,X);
% sh_dc = scnv(D,X_dc);
% sh_dnt = scnv(D,X_nt);

s_d = sh_d + sl;
% s_dc = sh_dc+sl_ref;
% s_dnt = sh_dnt + sl;

ps_nr = psnr(s_d,s_ref);
% ps_nr_dc = psnr(s_dc,s_ref);
% ps_nr_dnt = psnr(s_dnt,s_ref);

X_thre = X;
X_thre(X<0.15) = 0;
sh_thre = scnv(D,X_thre);
s_thre = sh_thre+sl;

figure;
imagesc(s_thre);colormap(gray);
title('threshold');
figure;
imagesc(s_d);colormap(gray);
title('original');


%plot various images
% figure;
% subplot(2,2,1);
% imagesc(sh_d); colormap(gray);
% title('hifreqdenoised');
% subplot(2,2,2);
% imagesc(sh); colormap(gray);
% title('noisyhighfreq');
% subplot(2,2,3);
% imagesc(s_d);colormap(gray);
% title(['final image psnr = ' num2str(ps_nr)]);
% subplot(2,2,4);
% imagesc(s_d - s_ref); colormap(gray);
% title('difference');
% 
% figure;
% subplot(2,2,1);
% imagesc(sh_dc); colormap(gray);
% title('Clean hifreqdenoised');
% subplot(2,2,2);
% imagesc(sh_ref); colormap(gray);
% title('Clean highfreq');
% subplot(2,2,3);
% imagesc(s_dc);colormap(gray);
% title(['Clean final image psnr = ' num2str(ps_nr_dc)]);
% subplot(2,2,4);
% imagesc(s_dc - s_ref); colormap(gray);
% title('Loss');
% 
% 
% figure;
% subplot(2,2,1);
% imagesc(sh_dnt); colormap(gray);
% title('Jnt hifreqdenoised');
% subplot(2,2,2);
% imagesc(sh); colormap(gray);
% title('Noisy highfreq');
% subplot(2,2,3);
% imagesc(s_dnt);colormap(gray);
% title(['Jnt final image psnr = ' num2str(ps_nr_dnt)]);
% subplot(2,2,4);
% imagesc(s_dc - s_ref); colormap(gray);
% title('Jnt Diff');
% 
% figure;
% subplot(2,2,1);
% imagesc(S_bm3dc);colormap(gray);
% title(['BM3D clean image psnr = ' num2str(ps_nr_bm3dc)]);
% subplot(2,2,2);
% imagesc(S_bm3dc - s_ref); colormap(gray);
% title('BM3D Clean Diff');
% subplot(2,2,3);
% imagesc(S_bm3d);colormap(gray);
% title(['BM3D denoised image psnr = ' num2str(ps_nr_bm3d)]);
% subplot(2,2,4);
% imagesc(S_bm3d - s_ref); colormap(gray);
% title('BM3D Diff');








