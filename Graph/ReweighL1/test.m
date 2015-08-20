% testing if reweighing L1 is feasible for removing white noise


%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load  Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load a saved noise
load('CacheData/stdnoise.mat');
sref = double(stdimage('lena.grey')) / 255;
sref = imresize(sref,.5);
[slref,shref] = lowpass(sref,5,15);
s = sref+r_noise;
[sl,sh] = lowpass(s,5,15);


% Load nice dictionary
load('CacheData/Dict_12x12.mat');




%%%%%%%%%%%%%%%%%%%%%%% set up parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = .05;
lambda = 0.3;
opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 100;
opt.rho = 10;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = Weight;


[Xcn,~] = cbpdn(D,sh,lambda,opt);
disp('cbpdn');

% opt.L1Weight = 1;
% [Xcnref,~] = cbpdn(D,shref,lambda,opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstructing conv
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

shreccn = scnv(D,Xcn);
sreccn = shreccn+sl;

% shreccnref = scnv(D,Xcnref);
% figure; imagesc(shreccnref); colormap(gray); colorbar;

p1 = psnr(sreccn,sref);
disp(['psnr from cbpdn: ',num2str(p1)]);


figure; 
imagesc(sreccn);
title(['cbpdn psnr = ',num2str(p1)]);
colormap(gray);
colorbar;


figure; 
imagesc(shreccn);
title(['cbpdn sh']);
colormap(gray);
colorbar;


figure; 
imagesc(sum3(Xcn));
title(['cbpdn coeff']);
colormap(gray);
colorbar;

