%script for testing white noise


% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('8x8x32');
num_dict = size(D,3);


%generate the noise
imsz = 256;
im_noise =single( 0.1*randn(imsz,imsz));
im_ref = single(zeros(imsz,imsz));

%setting up the parameters
lambda_all  = [0.1,0.15,0.2,0.25:0.03:0.52];
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

for i = 1:length(lambda_all)
    lambda = lambda_all(i);
    opt.rho = 100*lambda + 1;
    [X, optinf] = cbpdn(D, im_noise, lambda, opt);
    im_rec = convsum(D,X,1:1:num_dict);
    figure;
    imagesc(im_rec);
    colormap(gray); axis off;
    mse_rec(i) = immse(im_rec,im_ref);
    title(strcat('lambda = ',num2str(lambda),' mse = ',num2str(mse_rec(i))));
    saveas(gcf,strcat(sporco_path,'/Results/Denoising/wnoise/recon',num2str(i)),'png');
end

figure;
plot(lambda_all,mse_rec);
saveas(gcf,strcat(sporco_path,'/Results/Denoising/wnoise/curve'),'png');
