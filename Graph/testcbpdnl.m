
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load  Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load a saved noise
%load('noise_data.mat');
%snoise = s;
snoise = randn(256,256)*.1;
sref = single(stdimage('lena.grey')) / 255;
sref = imresize(sref,.5);
s = sref+snoise;
wsz = [20,20];
psz = [12,12];
neig = 30;
[sl,sh] = lowpass(s,5,15);
tau = 6;
[L,sh] = graphgen(sh,wsz,psz,neig,tau);
disp('graph generated');
sl = sl(1:size(sh,1),1:size(sh,2));

sref = sref(1:size(sh,1),1:size(sh,2));

[~,shref] = lowpass(sref,5,15);


% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');



%%%%%%%%%%%%%%%%%%%%%%% set up parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 1;
lambda = 0.25;
opt = {};
opt.Verbose = 0;
opt.MaxMainIter = 500;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
% load('wy.mat');
% y = 3*exp(-2*y);
%opt.L1Weight = repmat(y,[1,1,32]);
opt.L1Weight = 1;
%optimize
[Xnl,~]= cbpdn_L(D,sh,L,lambda,mu,opt);
disp('nl opt');

[Xcn,~] = cbpdn(D,sh,lambda,opt);
disp('cbpdn');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstructing conv
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');
shrecnl = scnv(D,Xnl);
shreccn = scnv(D,Xcn);

srecnl = shrecnl+sl;
sreccn = shreccn+sl;


p = psnr(srecnl,sref);
disp(['psnr from nonlocal: ',num2str(p)]);

p1 = psnr(sreccn,sref);
disp(['psnr from cbpdn: ',num2str(p1)]);

figure; 
imagesc(srecnl);
title(['nonlocal psnr = ',num2str(p)]);
colormap(gray);
colorbar;

figure; 
imagesc(sreccn);
title(['cbpdn psnr = ',num2str(p1)]);
colormap(gray);
colorbar;

figure; 
imagesc(shrecnl);
title(['nonlocal sh']);
colormap(gray);
colorbar;


figure; 
imagesc(shreccn);
title(['cbpdn sh']);
colormap(gray);
colorbar;

figure; 
imagesc(sum3(Xnl));
title(['nonlocal coeff']);
colormap(gray);
colorbar;

figure; 
imagesc(sum3(Xcn));
title(['cbpdn coeff']);
colormap(gray);
colorbar;




