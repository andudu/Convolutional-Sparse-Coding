%this script generates a noisy image and 

%batchname
batchname = 'lamplotsigma20';


% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('8x8x32');
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
lambda_all = 0.11:0.01:0.24;
mu_all = 0.5*lambda_all; %half the size of lambda
opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 300;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;

%vectors for record
ps_nr_all = zeros(size(lambda_all));
ps_nr_dnt_all = zeros(size(lambda_all));
ps_nr_bm3d_all = zeros(size(lambda_all));



%loop
disp('BM3D denoising');
[ps_nr_bm3d,S_bm3d] = BM3D(s_ref, s);
ps_nr_bm3d_all = ps_nr_bm3d*ones(size(lambda_all));    

for iter = 1:length(lambda_all)
    
    lambda = lambda_all(iter);
    mu = mu_all(iter);
    disp('jnt denoising');
    X_nt = cbpdnjnt(D,sh,lambda,mu,opt);
    disp('regular denoising');
    X = cbpdn(D,sh,lambda,opt);
    
    %computing
    sh_d = scnv(D,X);
    sh_dnt = scnv(D,X_nt);
    
    s_d = sh_d + sl;
    s_dnt = sh_dnt + sl;
    
    ps_nr = psnr(s_d,s_ref);
    ps_nr_dnt = psnr(s_dnt,s_ref);
    
    ps_nr_all(iter) = ps_nr;
    ps_nr_dnt_all(iter) = ps_nr_dnt;
    
    
    
    %plot various images
    figure;
    subplot(2,2,1);
    imagesc(sh_d); colormap(gray);
    title(['hifreqdenoised lambda = ',num2str(lambda)]);
    subplot(2,2,2);
    imagesc(sh); colormap(gray);
    title('noisyhighfreq');
    subplot(2,2,3);
    imagesc(s_d);colormap(gray);
    title(['final image psnr = ' num2str(ps_nr)]);
    subplot(2,2,4);
    imagesc(s_d - s_ref); colormap(gray);
    title('difference');
    saveas(gcf,[batchname,'_jntdenoiselam_',num2str(iter)],'fig');
    
    figure;
    subplot(2,2,1);
    imagesc(sh_dnt); colormap(gray);
    title(['Jnt hifreqd lambda = ',num2str(lambda)]);
    subplot(2,2,2);
    imagesc(sh); colormap(gray);
    title('Noisy highfreq');
    subplot(2,2,3);
    imagesc(s_dnt);colormap(gray);
    title(['Jnt final image psnr = ' num2str(ps_nr_dnt)]);
    subplot(2,2,4);
    imagesc(s_dc - s_ref); colormap(gray);
    title('Jnt Diff');
    saveas(gcf,[batchname,'_denoiselam_',num2str(iter)],'fig');
    
end

figure;
subplot(1,2,1);
imagesc(S_bm3d);colormap(gray);
title(['BM3D denoised image psnr = ' num2str(ps_nr_bm3d)]);
subplot(1,2,2);
imagesc(S_bm3d - s_ref); colormap(gray);
title('BM3D Diff');
saveas(gcf,[batchname,'_bm3d'],'fig');

figure;
plot(lambda_all,ps_nr_all,'b',lambda_all,ps_nr_dnt_all,'g',lambda_all,ps_nr_bm3d_all,'r');
legend('l1','jnt','bm3d');
saveas(gcf,[batchname,'_psnr'],'fig');













