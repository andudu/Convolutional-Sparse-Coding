%test script for super resolution( 2 step method)

% This sampling strategy actually causes a tiny shift  right and down!
% % Construct lowpass filtering and downsampling operators
% g = gauss2d([5 5], 0.75);
% smooth = @(x) imfilter(x, g, 'symmetric');
% dnsmpl = @(x) x(1:2:end, 1:2:end);
% dwnres = @(x) dnsmpl(smooth(x));

%Training on flicker images
flicker_ind = [1:1:20];
load('Flicker1_512_split.mat');
S0 = [];
k = 0;
for i = flicker_ind
    S0(:,:,(k)*4+1) = S(:,:,(i-1)*4+1);
    S0(:,:,(k)*4+2) = S(:,:,(i-1)*4+2);    
    S0(:,:,(k)*4+3) = S(:,:,(i-1)*4+3);
    S0(:,:,(k)*4+4) = S(:,:,(i-1)*4+4);
    k = k+1;
end
S0 = single(S0)/255;
image_num = length(flicker_ind);
clear S;

%Downsample image
S1 = zeros(128, 128, 4*image_num, 'single');
for k = 1:size(S0,3),
    S1(:,:,k) = imresize(S0(:,:,k),.5);
end

% Filter input images and compute highpass images
npd = 34;
fltlmbd = 5;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);
[Sl1, Sh1] = lowpass(S1, fltlmbd, npd);

% % Construct initial dictionary
numdict = 34;
D0 = zeros(12,12,numdict, 'single');
D0(4:9,4:9,:) = single(randn(6,6,numdict));
D0_bar = zeros(6,6,numdict, 'single');
D0_bar(2:5,2:5,:) = single(randn(4,4,numdict));

% load([sporco_path '/Data/ConvDict.mat']);
% dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
% f = dmap('8x8x32');
% D0 = f(:,:,1:12);

% Set up cbpdndliu parameters
lambda = 0.17;
opt = [];
opt.Verbose = 1;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(Sh,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.AutoDelta = 1;
opt.XRelaxParam = 1.5;
opt.DRelaxParam = 1.5;
opt.LinSolve = 'CG';

% Do dictionary learning
opt.MaxMainIter = 400;
[D_bar,X_bar, ~] = cbpdndliu(D0_bar, Sh1, lambda, opt);
X = PxT(X_bar);
[D,~] = ccmod(X,Sh,[12,12],opt);

tag = ['Flicker',num2str(image_num),'im',num2str(numdict),'dict','cold'];
% Display learned dictionary
o1.grey =1;
o1.unifscale =0;
square_plot(D,o1);
saveas(gcf,['Superres2stepResults/',tag,'_highresdict'],'fig');

% Display learned dictionary
o1.grey =1;
o1.unifscale =0;
square_plot(D_bar,o1);
saveas(gcf,['Superres2stepResults/',tag,'_lowresdict'],'fig');
o1.grey =0;
a = reshape(sum(abs(X_bar),3),size(X_bar,1),size(X_bar,2),size(X_bar,4));
square_plot(a,o1);

saveas(gcf,['Superres2stepResults/',tag,'_coeff'],'fig');

square_plot(Sh,o1);

square_plot(Sh1,o1);

save(['Superres2stepResults/',tag,'_dict.mat'],'D','D_bar');
save(['Superres2stepResults/',tag,'_coef.mat'],'X_bar');



% checknum = 13;
% figure;
% imagesc(convsum(D,X(:,:,:,checknum),1:1:13)); colormap(gray);
% title('high res rec');
% 
% figure;
% imagesc(Sh(:,:,checknum));colormap(gray);
% title('high res original');
% 
% figure;
% imagesc(convsum(D_bar,X_bar(:,:,:,checknum),1:1:13)); colormap(gray);
% title('low res rec');
% 
% figure;
% imagesc(Sh1(:,:,checknum)); colormap(gray);
% title('low res original');

