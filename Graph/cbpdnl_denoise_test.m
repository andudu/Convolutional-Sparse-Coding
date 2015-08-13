%script for testing and plotting cbpdn_L


%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load  Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load a saved noise
snoise = randn(256,256)*.1;
sref = double(stdimage('lena.grey')) / 255;
sref = imresize(sref,.5);
s = sref+snoise;
[sl,sh] = lowpass(s,5,15);


% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = double(dmap('12x12x36'));


%generate graph
optl = {};
optl.wsz = [60,60];
optl.psz = [12,12];
optl.neig = 50;
optl.Lformat = 'Sparse';
optl.Laplacian = 'n';
optl.Graph.tau = 2;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'Window';
optl.Graph.nsz = [7,7];
optl.Graph.k = [];
[L,sh] = laplacian_from_image(sh,optl);
disp('graph generated');

%cropping the image
sl = sl(1:size(sh,1),1:size(sh,2));
sref = sref(1:size(sh,1),1:size(sh,2));
[~,shref] = lowpass(sref,5,15);





%%%%%%%%%%%%%%%%%%%%%%% set up parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 2;
lambda = 0.24;
opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 50;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;
opt.Ysolver = 'fista';


%optimize
[Xnl,~]= cbpdn_L(D,sh,L,lambda,mu,opt);
disp('nl opt');

opt = rmfield(opt,'Ysolver');
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




