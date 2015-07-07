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
lambda = 0.1;
mu_all = [0,0.02:0.005:0.055];
ps_nr_dnt_all = zeros(size(mu_all));
opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 300;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;


batchname = 'firstbatch';

for iter = 1:length(mu_all)
mu = mu_all(iter);
disp('jnt denoising');
X_nt = cbpdnjnt(D,sh,lambda,mu,opt);
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
saveas(gcf,[batchname,'_jntdenoise_',num2str(iter)],'fig');
end
figure;
plot(mu_all,ps_nr_dnt_all);
saveas(gcf,['psnr_mu_',batchname],'fig');
save([batchname,'.mat'],'mu_all','ps_nr_dnt_all');










