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




opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 400;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
lambda = 0.29; % best lambda from parameter search
mu_grim = []; %best mu's
mu_grcoef = [];

w_all = 0.5:0.05:1.5;
snr_all_grim = zeros(size(w_all)); %record both gr im and gr coeff
snr_all_grcoef = zeros(size(w_all));
psnr_all_grim = zeros(size(w_all));
psnr_all_grcoef = zeros(size(w_all));
X_max_grim = [];
mu_max_grim = 0;
mu_max_grcoef = 0;
psnr_max_grim = 0;
psnr_max_grcoef = 0;


disp('Running L1 weight search');
for iter = 1:length(w_all)
%setting up parameters
w = w_all(iter);
opt.L1Weight =reshape([ones(1,num_dict-1),w],1,1,num_dict);  %simple L1 first
disp(['iter: ',num2str(iter)]);
mu = mu_all(iter);
[X_grim, ~] = cbpdngr_new(D, s, lambda, mu_grim, opt);
[X_grcoef, ~] = cbpdngr(D, s, lambda, mu_grcoef, opt);
DX_grim = scnv(D,X_grim);
DX_grcoef = scnv(D,X_grcoef);


snr_all_grim(iter) = snr(s_ref,DX_grim);
psnr_all_grim(iter) = psnr(DX_grim,s_ref);
snr_all_grcoef(iter) = snr(s_ref,DX_grcoef);
psnr_all_grcoef(iter) = psnr(DX_grcoef,s_ref);



if(psnr_all_grim(iter)> psnr_max_grim)
    psnr_max_grim = psnr_all_grim(iter);
    w_max_grim = w_all(iter);
    X_max_grim = X_grim;
end


if(psnr_all_grcoef(iter)> psnr_max_grcoef)
    psnr_max_grcoef = psnr_all_grcoef(iter);
    w_max_grcoef = w_all(iter);
    X_max_grcoef = X_grcoef;
end


end

%save the info for the max point
save(strcat(sporco_path,'/Results/Denoising/Raw/mu_search_nofilter_',...
    imagename,'.mat'),'X_max_grim','X_max_grcoef','psnr_all_grim',...
    'psnr_all_grcoef','psnr_max_grim','psnr_max_grcoef',...
    'snr_all_grim','snr_all_grcoef','mu_grim','mu_grcoef'...
    ,'w_all','w_max_grim','w_max_grcoef','lambda');


%plot and save the best reconstruction(grim)
DX_max_grim = scnv(D,X_max_grim);
diff_max_grim = DX_max_grim - s_ref;

figure;
subplot(1,2,1);
imagesc(DX_max_grim); colormap (gray); axis off;
title(strcat('Recon:  wgrim = ',num2str(w_max_grim),'  psnr = ',num2str(psnr_max_grim)));
subplot(1,2,2);
imagesc(diff_max_grim); colormap (gray); axis off;
title('Difference');
saveas(gcf,strcat(sporco_path,'/Results/Denoising/best_w_search_grim_',imagename),'png');

%plot and save the best reconstruction(grcoef)
DX_max_grcoef = scnv(D,X_max_grcoef);
diff_max_grcoef = DX_max_grcoef - s_ref;

figure;
subplot(1,2,1);
imagesc(DX_max_grcoef); colormap (gray); axis off;
title(strcat('Recon:  wgrcoef = ',num2str(mu_max_grcoef),'  psnr = ',num2str(psnr_max_grcoef)));
subplot(1,2,2);
imagesc(diff_max_grcoef); colormap (gray); axis off;
title('Difference');
saveas(gcf,strcat(sporco_path,'/Results/Denoising/best_w_search_grcoef_',imagename),'png');



