%this script compares sparse coding between applying the filter 
%and ways to avoid it. 

% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
num_dict = size(D,3);

imagename = 'lena';

% Load test image
s = single(rgbtogrey(stdimage(imagename)))/255; %use barbara or lena here
s_ref = s;
%add the noise (10 percent peak signal)
sigma = .1;
s = s+ sigma*randn(size(s));

if isempty(s),
  error('Data required for demo scripts has not been installed.');
end

scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

% Highpass filter test image
npd = 16;
fltlmbd = 5;
[sl, sh] = lowpass(s, fltlmbd, npd);

lambda_all = 0.01:0.02:0.55;
snr_all = zeros(size(lambda_all));
psnr_all = zeros(size(lambda_all));
X_max = [];
optinf_max = {};
lambda_max = 0;
psnr_max = 0;


for iter = 1:length(lambda_all)
%setting up parameters
disp(['iter: ',num2str(iter)]);
lambda = lambda_all(iter);
opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 400;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;
[X1, optinf1] = cbpdn(D, sh, lambda, opt);
%save(strcat(sporco_path,'/Results/benchmark_',imagename,'.mat'),'X1','optinf1');
H1 = scnv(D,X1);
DX1 = H1+sl;

snr_all(iter) = snr(s_ref,DX1);
psnr_all(iter) = psnr(DX1,s_ref);

if(psnr_all(iter)> psnr_max)
    psnr_max = psnr_all(iter);
    lambda_max = lambda_all(iter);
    X_max = X1;
    optinf_max = optinf1;
end



end

%save the info for the max point
save(strcat(sporco_path,'/Results/Denoising/Raw/lambda_search_filter_',...
    imagename,'.mat'),'X_max','optinf_max','psnr_all','psnr_max'...
    ,'lambda_all','lambda_max','snr_all');

%plot and save the best reconstruction
H_max = scnv(D,X_max);
DX_max = H_max+sl;
diff_max = DX_max - s_ref;

figure;
subplot(1,2,1);
imagesc(DX_max); colormap (gray); axis off;
title(strcat('Recon:  lambda_{max} = ',num2str(lambda_max),'  snr = ',num2str(psnr_max)));
subplot(1,2,2);
imagesc(diff_max); colormap (gray); axis off;
title('Difference');
saveas(gcf,strcat(sporco_path,'/Results/Denoising/best_lambda_search_',imagename),'png');








