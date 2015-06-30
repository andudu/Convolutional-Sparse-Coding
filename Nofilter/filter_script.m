%this script compares sparse coding between applying the filter 
%and ways to avoid it. 

% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('8x8x32');
num_dict = size(D,3);

imagename = 'lena';

% Load test image
s = single(rgbtogrey(stdimage(imagename)))/255; %use barbara here
if isempty(s),
  error('Data required for demo scripts has not been installed.');
end

scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

% Highpass filter test image
npd = 16;
fltlmbd = 5;
[sl, sh] = lowpass(s, fltlmbd, npd);



mu_all = 0.01:0.01:0.06;


%setting up parameters
lambda = 0.01;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 500;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;
[X1, optinf1] = cbpdn(D, sh, lambda, opt);
%save(strcat(sporco_path,'/Results/benchmark_',imagename,'.mat'),'X1','optinf1');
H1 = scnv(D,X1);
DX1 = H1+sl;
res1 = snr(s,DX1);
diff_1 = s-DX1;

for iter = 1:length(mu_all)

% Compute representation using cbpdngr with no pass but high regularizer
mu = mu_all(iter);
opt.L1Weight = reshape([ones(1,num_dict-1),.8],1,1,num_dict); %try few sparsity .9 for new_grd
opt.GrdWeight = reshape([zeros(1,num_dict-1),1],1,1,num_dict); %5 on coeff 1 for final
[X2, optinf2] = cbpdngr_new(D, s, lambda, mu, opt);

save(strcat(sporco_path,'/Results/Raw/new_grdx_',imagename,num2str(iter),'.mat'),'X2','optinf2');

%plot reconstruction 
DX2 = scnv(D,X2);
res2 = snr(s,DX2);
figure;
subplot(1,2,1);
imagesc(DX1); colormap (gray); axis off;
title(strcat('lambda = ',num2str(lambda),'  snr = ',num2str(res1)));
subplot(1,2,2);
imagesc(DX2); colormap (gray); axis off;
title(strcat('mu = ',num2str(mu),'  snr = ',num2str(res2)));
saveas(gcf,strcat(sporco_path,'/Results/newgrd_full_rec_',imagename,num2str(iter)),'png');


%plot high freq rec
diff_2 = convsum(D,X2,1:1:num_dict-1);
figure;
subplot(1,2,1);
imagesc(H1); colormap (gray); axis off;
title(strcat('lambda = ',num2str(lambda)));
subplot(1,2,2);
imagesc(diff_2); colormap (gray); axis off;
title(strcat('mu = ',num2str(mu)));
saveas(gcf,strcat(sporco_path,'/Results/newgrd_High_rec_',imagename,num2str(iter)),'png');


%plot low freq 
L2 = convsum(D,X2,num_dict);
figure;
subplot(1,2,1);
imagesc(sl); colormap (gray); axis off;
title(strcat('lambda = ',num2str(lambda)));
subplot(1,2,2);
imagesc(L2); colormap (gray); axis off;
title(strcat('mu = ',num2str(mu)));
saveas(gcf,strcat(sporco_path,'/Results/newgrd_Low_rec_',imagename,num2str(iter)),'png');

%plot difference
diff_2 = s-DX2;
figure;
subplot(1,2,1);
imagesc(diff_1); colormap (gray); axis off;
title(strcat('lambda = ',num2str(lambda)));
subplot(1,2,2);
imagesc(diff_2);  axis off;
title(strcat('mu = ',num2str(mu)));
saveas(gcf,strcat(sporco_path,'/Results/newgrd_Diff_rec_',imagename,num2str(iter)),'png');


end






