
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load  Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load a saved noise
load('noise_data.mat');
snoise = s;
% snoise = randn(256,256)*.1;
s_ref = single(stdimage('lena.grey')) / 255;
s_ref = imresize(s_ref,.5);
s = s_ref+snoise;
wsz = [60,60];
psz = [12,12];
neig = 30;
[sl,sh] = lowpass(s,5,15);
[L,sh] = graphgen(sh,wsz,psz,neig);
disp('graph generated');
sl = sl(1:size(sh,1),1:size(sh,2));


% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');



%%%%%%%%%%%%%%%%%%%%%%% set up parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = .6;
opt = {};
opt.Verbose = 0;
opt.MaxMainIter = 500;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;

%optimize
lambda = 0.24;
[Xnl,~]= cbpdn_L(D,sh,L,lambda,mu,opt);
disp('nl opt');
lambda = 0.27;
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

sref = sref(1:size(sh,1),1:size(sh,2));

p = psnr(srecnl,sref);
disp(['psnr from nonlocal: ',num2str(p)]);

p1 = psnr(sreccn,sref);
disp(['psnr from cbpdn: ',num2str(p1)]);

figure; 
imagesc(srecnl);
title(['nonlocal psnr = ',num2str(p)]);
colorbar;

figure; 
imagesc(sreccn);
title(['cbpdn psnr = ',num2str(p1)]);
colorbar;

figure; 
imagesc(shrecnl);
title(['nonlocal sh']);
colorbar;


figure; 
imagesc(shreccn);
title(['cbpdn sh']);
colorbar;

figure; 
imagesc(sum3(Xnl));
title(['nonlocal coeff']);
colorbar;

figure; 
imagesc(sum3(Xcn));
title(['cbpdn coeff']);
colorbar;




