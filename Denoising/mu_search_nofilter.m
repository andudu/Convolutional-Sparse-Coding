%this script compares sparse coding between applying the filter 
%and ways to avoid it. 

% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('8x8x32');
num_dict = size(D,3);

tag = 'gauss' ;%pointsource


%resize
fsz = size(D,1);
gsz = 40;
npad = gsz - fsz;
D = padarray(D,[npad,npad],'post');




% %add this gaussian
D(:,:,num_dict) = gauss2d(gsz,10);


% %add point source
% pixsz = 1;
% D_point = ones(pixsz,pixsz);
% D_point = padarray(D_point,[gsz-pixsz,gsz-pixsz],'post');
% D(:,:,num_dict) = D_point;

%normalize
D(:,:,num_dict) = D(:,:,num_dict)/norm(vec(D(:,:,num_dict)));





imagename = 'lena';

% Load test image
s = single(rgbtogrey(stdimage(imagename)))/255; %use barbara or lena here
s = imresize(s, 0.5);
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


lambda = .2; % best lambda from parameter search

opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 400;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight =1;  %simple L1 first
opt.AutoRhoPeriod = 10;
opt.RelaxParam = 1;
disp('Best Reconstruction from L1 with filter');
% [X_ref, ~] = cbpdn(D, sh, lambda, opt);
% H_ref = scnv(D,X_ref);
% DX_ref = H_ref+sl;
% psnr_ref = psnr(DX_ref,s_ref); %psnr for without gradients
opt.GrdWeight = reshape([zeros(1,num_dict-1),1],1,1,num_dict);







mu_all = 5;

w_all =5*10e-5;



snr_all_grim = zeros(length(mu_all),length(w_all)); %record both gr im and gr coeff
snr_all_grcoef = zeros(length(mu_all),length(w_all));
psnr_all_grim = zeros(length(mu_all),length(w_all));
psnr_all_grcoef = zeros(length(mu_all),length(w_all));
X_max_grim = [];
mu_max_grim = 0;
w_max_grim = 0;
mu_max_grcoef = 0;
psnr_max_grim = 0;
psnr_max_grcoef = 0;

disp('Running mu w search');
for iter = 1:length(mu_all)
    mu = mu_all(iter);
    for iter2 = 1:length(w_all)
        %setting up parameters
        w = w_all(iter2);
        opt.L1Weight = reshape([ones(1,num_dict-1),w],1,1,num_dict);
        disp(['iter: ',num2str(length(mu_all)*(iter-1)+iter2)]);
        %         %%%%%%%%%%%%%%%%%%%%%
        %         [X_grcoef, ~] = cbpdngr(D, s, lambda, mu, opt);
        %         DX_grcoef = scnv(D,X_grcoef);
        %
        %         snr_all_grcoef(iter,iter2) = snr(s_ref,DX_grcoef);
        %         psnr_all_grcoef(iter,iter2) = psnr(DX_grcoef,s_ref);
        %
        %
        %         if(psnr_all_grcoef(iter,iter2)> psnr_max_grcoef)
        %             psnr_max_grcoef = psnr_all_grcoef(iter,iter2);
        %             mu_max_grcoef = mu_all(iter,iter2);
        %             X_max_grcoef = X_grcoef;
        %         end
        %         figure;
        %         subplot(1,3,1);
        %         imagesc(convsum(D,X_grcoef,1:1:num_dict)); colormap (gray); axis off;
        %         title(strcat('GrcoefHigh: mu = ',num2str(mu),' w = ',num2str(w) ));
        %         subplot(1,3,2);
        %         imagesc(convsum(D,X_grcoef,num_dict)); colormap (gray); axis off;
        %         title(strcat('GrcoefLow: mu = ',num2str(mu),' w = ',num2str(w) ));
        %         saveas(gcf,strcat(sporco_path,'/Results/Denoising/grcoef_',num2str(iter),num2str(iter2),imagename),'png');
        %         subplot(1,3,3);
        %         imagesc(DX_grcoef); colormap (gray); axis off;
        %         title(strcat('GrcoefRec: mu = ',num2str(mu),' w = ',num2str(w),'  psnr = ',num2str(psnr_all_grcoef(iter,iter2))));
        
        %%%%%%%%%%%%%%%%%%%%%
       
         [X_grim, ~] = cbpdngr_new(D, s, lambda, mu, opt);
%         DX_grim = scnv(D,X_grim);
        [X_grim, Z,~] = cbpdnlc(D, s, lambda, mu, opt);
        snr_all_grim(iter,iter2) = snr(s_ref,DX_grim);
        psnr_all_grim(iter,iter2) = psnr(DX_grim,s_ref);
        
        if(psnr_all_grim(iter,iter2)> psnr_max_grim)
            psnr_max_grim = psnr_all_grim(iter,iter2);
            mu_max_grim = mu_all(iter);
            w_max_grim = w_all(iter2);
            X_max_grim = X_grim;
        end
        figure;
        subplot(1,3,1);
        imagesc(convsum(D,X_grim,1:1:num_dict-1)); colormap (gray); axis off;
        title(strcat('GrimHigh: mu = ',num2str(mu),' w = ',num2str(w) ));
        subplot(1,3,2);
        imagesc(convsum(D,X_grim,num_dict)); colormap (gray); axis off;
        title(strcat('GrimLow: mu = ',num2str(mu),' w = ',num2str(w) ));
        subplot(1,3,3);
        imagesc(DX_grim); colormap (gray); axis off;
        title(strcat('GrimRec: mu = ',num2str(mu),' w = ',num2str(w),'  psnr = ',num2str(psnr_all_grim(iter,iter2)) ));
        saveas(gcf,strcat(sporco_path,'/Results/Denoising/grim_',num2str(iter),num2str(iter2),imagename,tag),'png');
        %
    end
end

%save the info for the max point
save(strcat(sporco_path,'/Results/Denoising/Raw/mu_search_nofilter_',...
    imagename,tag,'.mat'),'X_max_grim','psnr_all_grim','psnr_max_grim',...
    'snr_all_grim'...
    ,'mu_all','mu_max_grim','lambda');
% save(strcat(sporco_path,'/Results/Denoising/Raw/mu_search_nofilter_',...
%     imagename,'.mat'),'X_max_grim','X_max_grcoef','psnr_all_grim',...
%     'psnr_all_grcoef','psnr_max_grim','psnr_max_grcoef',...
%     'snr_all_grim','snr_all_grcoef'...
%     ,'mu_all','mu_max_grim','mu_max_grcoef','lambda');


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
saveas(gcf,strcat(sporco_path,'/Results/Denoising/best_mu_search_grim_',imagename,tag),'png');

% %plot and save the best reconstruction(grcoef)
% DX_max_grcoef = scnv(D,X_max_grcoef);
% diff_max_grcoef = DX_max_grcoef - s_ref;
% 
% figure;
% subplot(1,2,1);
% imagesc(DX_max_grcoef); colormap (gray); axis off;
% title(strcat('Recon:  mugrcoef = ',num2str(mu_max_grcoef),'  psnr = ',num2str(psnr_max_grcoef)));
% subplot(1,2,2);
% imagesc(diff_max_grcoef); colormap (gray); axis off;
% title('Difference');
% saveas(gcf,strcat(sporco_path,'/Results/Denoising/best_mu_search_grcoef_',imagename),'png');





