%this script generates a noisy image and 

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
%add the noise (10 percent peak signal)
sigma = 15/255;
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
lambda = 0.1335; %best lambda
mu_all = [0,0.001:0.001:0.009,0.01:0.01:0.09,0.1:0.1:0.5];
ps_nr_dnt_all = zeros(size(mu_all));
opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 300;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;


%compute the best psnr for regular cbpdn
disp('Regular CBPDN');
[X_d, ~] = cbpdn(D, sh, lambda, opt);
sh_d = scnv(D,X_d);
s_d = sh_d + sl;
ps_nr_d = psnr(s_d,s_ref);

%plot various images
figure;
subplot(2,2,1);
imagesc(sh_d); colormap(gray);
title('hifreqdenoised');
subplot(2,2,2);
imagesc(sh); colormap(gray);
title('noisyhighfreq');
subplot(2,2,3);
imagesc(s_d);colormap(gray);
title(['final image psnr = ' num2str(ps_nr_d)]);
subplot(2,2,4);
imagesc(s_d - s_ref); colormap(gray);
title('difference');
saveas(gcf,[sporco_path,'/Results/Denoising/GroupNorm/',batchname,'_Regdenoise'],'fig');

batchname = 'MuSearchSigma15New';



%loop through the joint norm cbpdn
for iter = 1:length(mu_all)
mu = mu_all(iter);
disp('jnt denoising');
X_nt = cbpdnjnt_new(D,sh,lambda,mu,opt);
%computing 
sh_dnt = scnv(D,X_nt);
s_dnt = sh_dnt + sl;
ps_nr_dnt = psnr(s_dnt,s_ref);
ps_nr_dnt_all(iter) = ps_nr_dnt;
figure;
subplot(2,2,1);
imagesc(sh_dnt); colormap(gray);
title('Jnt hifreqdenoised');
subplot(2,2,2);
imagesc(sh_dnt); colormap(gray);
title('Jnt highfreq');
subplot(2,2,3);
imagesc(s_dnt);colormap(gray);
title(['Jnt final image psnr = ' num2str(ps_nr_dnt)]);
subplot(2,2,4);
imagesc(s_dc - s_ref); colormap(gray);
title(['Jnt Diff mu = ',num2str(mu)]);
saveas(gcf,[sporco_path,'/Results/Denoising/GroupNorm/',batchname,'_jntdenoise_',num2str(iter)],'fig');
end
figure;
plot(log(mu_all),ps_nr_dnt_all,'b');
hold on;
plot(log(mu_all), ps_nr_d*ones(size(mu_all)),'r');
legend('joint','regular');
title('sigma = 15/255');

saveas(gcf,[sporco_path,'/Results/Denoising/GroupNorm/','psnr_mu_',batchname],'fig');
save([sporco_path,'/Results/Denoising/GroupNorm/',batchname,'.mat'],'mu_all','ps_nr_dnt_all','sigma');










