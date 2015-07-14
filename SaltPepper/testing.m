% Script demonstrating usage of the varsplit_test function.


% Training images
s = imresize(single(stdimage('lena.grey'))/255,.5) ;


% Add Salt and Pepper
% Sn = S0;
% Sn(250:260,250:260) = 1;
% s = imnoise(s, 'salt & pepper');




% Highpass filter test image

npd = 16;
fltlmbd = 5;
[sl, sh] = lowpass(s, fltlmbd, npd);
shn = sh;
opt = [];
I = randi(255,2,ceil(256*256*.08));
opt.L1Weight = ones(256,256,33);
for i = 1:size(I,2)
    opt.L1Weight(I(1,i),I(2,i),33) = 0;
    shn(I(1,i),I(2,i)) = 1;
end


% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('8x8x32');
delta = zeros(8,8);
delta(1,1) = 1;
D(:,:,end+1) = delta;


lambda = .02;
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.AutoRhoPeriod = 10;


[X,~] = cbpdn(D,shn,lambda,opt);
% opt.L1Weight = 1;
% [X_clean,~] = cbpdn(D(:,:,1:1:32),Sh,lambda,opt);



% figure;
% title('Lowpass');
% imagesc(Sln);

Sh_rec = convsum(D,X,1:1:32);

figure;
imagesc(shn);
title('noisy');

figure;
imagesc(Sh_rec);
title('Reconstructed High Freq');

figure;
imagesc(sh);
title('Original');

% figure;
% p = psnr(Sh_clean+Sl,S0);
% imagesc(Sl+Sh_clean);
% title(['psnr of clean recreation = ' num2str(p)]);




% temp = zeros(size(X,1),size(X,2),32);
% for i = 1:32
%     temp(:,:,i) = convsum(D,X,i);
%     figure;
%     imagesc(temp(:,:,i));
% end
     









