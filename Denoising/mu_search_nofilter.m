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



opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 400;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight =1;  %simple L1 first
disp('Best Reconstruction from L1 with filter');
lambda = 0.29; % best lambda from parameter search
[X_ref, ~] = cbpdn(D, sh, lambda, opt);
H_ref = scnv(D,X_ref);
DX_ref = H_ref+sl;
psnr_ref = psnr(DX_ref,s_ref); %psnr for without gradients
opt.GrdWeight = reshape([zeros(1,num_dict-1),1],1,1,num_dict);


mu_all = 0.01:0.02:0.53;

snr_all_grim = zeros(size(mu_all)); %record both gr im and gr coeff
snr_all_grcoef = zeros(size(mu_all));
psnr_all_grim = zeros(size(mu_all));
psnr_all_grcoef = zeros(size(mu_all));
X_max_grim = [];
mu_max_grim = 0;
mu_max_grcoef = 0;
psnr_max_grim = 0;
psnr_max_grcoef = 0;

disp('Running mu search');
for iter = 1:length(mu_all)
%setting up parameters
disp(['iter: ',num2str(iter)]);
mu = mu_all(iter);
[X_grim, ~] = cbpdngr_new(D, s, lambda, mu, opt);
[X_grcoef, ~] = cbpdngr(D, s, lambda, mu, opt);
DX_grim = scnv(D,X_grim);
DX_grcoef = scnv(D,X_grcoef);


snr_all_grim(iter) = snr(s_ref,DX_grim);
psnr_all_grim(iter) = psnr(DX_grim,s_ref);
snr_all_grcoef(iter) = snr(s_ref,DX_grcoef);
psnr_all_grcoef(iter) = psnr(DX_grcoef,s_ref);



if(psnr_all_grim(iter)> psnr_max_grim)
    psnr_max_grim = psnr_all_grim(iter);
    mu_max_grim = mu_all(iter);
    X_max_grim = X_grim;
end


if(psnr_all_grcoef(iter)> psnr_max_grcoef)
    psnr_max_grcoef = psnr_all_grcoef(iter);
    mu_max_grcoef = mu_all(iter);
    X_max_grcoef = X_grcoef;
end


end

%save the info for the max point
save(strcat(sporco_path,'/Results/Denoising/Raw/mu_search_nofilter_',...
    imagename,'.mat'),'X_max_grim','X_max_grcoef','psnr_all_grim',...
    'psnr_all_grcoef','psnr_max_grim','psnr_max_grcoef','psnr_ref',...
    'snr_all_grim','snr_all_grcoef','X_ref'...
    ,'mu_all','mu_max_grim','mu_max_grcoef','lambda');


%plot and save the best reconstruction(grim)
DX_max_grim = scnv(D,X_max_grim);
diff_max_grim = DX_max_grim - s_ref;

figure;
subplot(1,2,1);
imagesc(DX_max_grim); colormap (gray); axis off;
title(strcat('Recon:  mugrim = ',num2str(mu_max_grim),'  psnr = ',num2str(psnr_max_grim)));
subplot(1,2,2);
imagesc(diff_max_grim); colormap (gray); axis off;
title('Difference');
saveas(gcf,strcat(sporco_path,'/Results/Denoising/best_mu_search_grim_',imagename),'png');

%plot and save the best reconstruction(grcoef)
DX_max_grcoef = scnv(D,X_max_grcoef);
diff_max_grcoef = DX_max_grcoef - s_ref;

figure;
subplot(1,2,1);
imagesc(DX_max_grcoef); colormap (gray); axis off;
title(strcat('Recon:  mugrcoef = ',num2str(mu_max_grcoef),'  psnr = ',num2str(psnr_max_grcoef)));
subplot(1,2,2);
imagesc(diff_max_grcoef); colormap (gray); axis off;
title('Difference');
saveas(gcf,strcat(sporco_path,'/Results/Denoising/best_mu_search_grcoef_',imagename),'png');





