%script for testing white noise


% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('8x8x32');
D2 = dmap('12x12x36');
num_dict = size(D,3);


%generate the noise
imsz = 256;
im_ref = single(zeros(imsz,imsz));
im_ref(126:128,:) = 1; %add a line
im_noise =im_ref+single( 0.07*randn(imsz,imsz));


%setting up the parameters
% lambda_all  = [0.1,0.15,0.2,0.25:0.03:0.52];
lambda_all = 0.5;
opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 300;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.AutoRhoPeriod = 10;
opt.RelaxParam = 1;


%record the variables
mse_rec = zeros(size(lambda_all));


%original image
figure;
imagesc(im_noise);
colormap(gray);
axis off;
title('original image');

for i = 1:length(lambda_all)
    lambda = lambda_all(i);
    opt.rho = 100*lambda + 1;
    [X1, optinf] = cbpdn(D, im_noise, lambda, opt);
    im_rec = convsum(D,X,1:1:num_dict);
    figure;
    imagesc(im_rec);
    colormap(gray); axis off;
    mse_rec(i) = immse(im_rec,im_ref);
    title(strcat('SmallD lambda = ',num2str(lambda),' mse = ',num2str(mse_rec(i))));
    saveas(gcf,strcat(sporco_path,'/Results/Denoising/wnoise/recon1',num2str(i)),'png');
    figure;
    imagesc(im_rec-im_ref);
    colormap(gray);
    
    
    [X2, optinf] = cbpdn(D2, im_noise, lambda, opt);
    im_rec = convsum(D2,X,1:1:num_dict);
    figure;
    imagesc(im_rec);
    colormap(gray); axis off;
    mse_rec(i) = immse(im_rec,im_ref);
    title(strcat('LargeD lambda = ',num2str(lambda),' mse = ',num2str(mse_rec(i))));
    saveas(gcf,strcat(sporco_path,'/Results/Denoising/wnoise/recon2',num2str(i)),'png');  
    figure;
    imagesc(im_rec-im_ref);
    colormap(gray);
    
end

% figure;
% plot(lambda_all,mse_rec);
% saveas(gcf,strcat(sporco_path,'/Results/Denoising/wnoise/curve'),'png');
