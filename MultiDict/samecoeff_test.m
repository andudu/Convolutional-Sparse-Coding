% Script demonstrating usage of the varsplit_test function.

% 
% % Training images
% S0 = zeros(512, 512, 5, 'single');
% S0(:,:,1) = single(stdimage('lena.grey')) / 255;
% S0(:,:,2) = single(stdimage('barbara.grey')) / 255;
% S0(:,:,3) = single(stdimage('kiel.grey')) / 255;
% S0(:,:,4) = single(rgb2gray(stdimage('mandrill'))) / 255;
% tmp = single(stdimage('man.grey')) / 255;
% S0(:,:,5) = tmp(101:612, 101:612);
% 
% 
% %Reduce images size to speed up demo script
% tmp = zeros(256, 256, 5, 'single');
% for k = 1:size(S0,3),
%   tmp(:,:,k) = imresize(S0(:,:,k), 0.5);
% end
% S1 = tmp;

% %Reduce images size to speed up demo script
% tmp = zeros(128, 128, 5, 'single');
% for k = 1:size(S0,3),
%   tmp(:,:,k) = imresize(S0(:,:,k), 0.5);
% end
% S1 = tmp;


%Training on flicker images
flicker_ind = [4,5,7,10,14,16,17,19,23,30,43,46,49,30,35,36];
load('FlickrCC_512_512.mat');
%flicker_ind = flicker_ind(1:5);
S0 = single(S(:,:,flicker_ind))/255;
image_num = length(flicker_ind);
clear S;


% %Reduce images size to speed up demo script
% tmp = zeros(256, 256, image_num, 'single');
% for k = 1:size(S0,3),
%   tmp(:,:,k) = imresize(S0(:,:,k), 0.5);
% end
% S1 = tmp;

%Reduce images size to speed up demo script
tmp = zeros(512, 512, image_num, 'single');
for k = 1:size(S0,3),
  a = imresize(S0(:,:,k), 0.5);
  tmp(:,:,k) = imresize(a, 2);
end
S1 = tmp;

% Filter input images and compute highpass images
npd = 16;
fltlmbd = 6;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);
[Sl1, Sh1] = lowpass(S1, fltlmbd, npd);
%[Sl2,Sh2] = lowpass(S2,fltlmbd,npd);

% % Construct initial dictionary
numdict = 32;
D0 = zeros(12,12,numdict, 'single');
D0(4:9,4:9,:) = single(randn(6,6,numdict));
% load([sporco_path '/Data/ConvDict.mat']);
% dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
% f = dmap('8x8x32');
% D0 = f(:,:,1:12);

% Set up cbpdndliu parameters
lambda = 0.2;
opt = [];
opt.Verbose = 1;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(Sh,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.AutoDelta = 1;

%just add a seperate stepsize for all _bar variable
opt.AutoRho_bar = 1;
opt.AutoRhoPeriod_bar = 10;
opt.AutoSigma_bar = 1;
opt.AutoSigmaPeriod_bar = 10;
opt.AutoDelta = 1;
opt.AutoDelta = 10;
opt.XRelaxParam = 1.4;
opt.DRelaxParam = 1.4;

% Do dictionary learning
opt.MaxMainIter = 250;
[D, D_bar,X, optinf] = samecoeff_multilearn(D0, Sh, Sh1, lambda, opt);
reso = 512;
tag = ['Flicker',num2str(image_num),'im',num2str(numdict),'dict',num2str(reso),'cold'];
% Display learned dictionary
o1.grey =1;
o1.unifscale =0;
square_plot(D,o1);
saveas(gcf,['SameCoeffResults/',tag,'_highresdict'],'fig');

% Display learned dictionary
o1.grey =1;
o1.unifscale =0;
square_plot(D_bar,o1);
saveas(gcf,['SameCoeffResults/',tag,'_lowresdict'],'fig');
o1.grey =0;
a = reshape(sum(abs(X),3),size(X,1),size(X,2),size(X,4));
square_plot(a,o1);

saveas(gcf,['SameCoeffResults/',tag,'_coeff'],'fig');

save(['SameCoeffResults/',tag,'_dict.mat'],'D','D_bar','X');

% num = 2;
% a = reshape(sum(abs(X(:,:,num,:)),3),size(X,1),size(X,2),size(X,4));
% square_plot(a,o1);
% a = reshape(sum(abs(X_bar(:,:,num,:)),3),size(X_bar,1),size(X_bar,2),size(X_bar,4));
% square_plot(a,o1);


% for i = 1:5
%     temp = convsum(D,X(:,:,:,i),1:14);
%     figure;
%     imagesc(temp);
%     temp = convsum(D_bar,X(:,:,:,i),1:14);
%     figure;
%     imagesc(temp);
%     figure;
%     imagesc(Sh(:,:,i));
%     
% end


% Plot functional value evolution
figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');