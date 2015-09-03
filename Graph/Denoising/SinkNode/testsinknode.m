% script for testing sink node graph for Laplacian Sparse Coding

% tested using lena patch

% small scale denoise test 

%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting options for Graph Laplacian
optl = {};
optl.wsz = [60,60];
optl.psz = [12,12];
optl.neig = 40;
optl.Lformat = 'Full';
optl.Laplacian = 'u';
optl.Graph.tau = 1;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'Window';
optl.Graph.nsz = [9,9];
optl.Graph.k = [];
optl.Graph.WeightMatrix = 1;

imsz = optl.wsz;
psz = optl.psz;
stpsz = [1,1];
load('stdnoise.mat');
n = r_noise(1:imsz(1),1:imsz(2));
Dversion = 'lenapatch';
%%%%%%%%%%%%%%%%%%   Generating Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(Dversion, 'lenapatch'),
    %patch from lena
    s = double(stdimage('lena.grey')) / 255;
    s = imresize(s,.5);
    s = s+.1*randn(size(s));
    [sl,sh] = lowpass(s,4,15);
    sh = sh(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
    sl = sl(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
end



    


%%%%%%%%%%%%%%%%%   Building the Graph Weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L,scrop] = laplacian_from_image(sh,optl);
Ltemp = L{1}.M;
Ltemp(1:size(Ltemp,1)+1:end) = 0;
d = sum(Ltemp,1); %degree
Ltemp(1:size(Ltemp,1)+1:end) = -d;
Ltemp = Ltemp/mean(d);  %naive normalization

d = d./(.2*max(d));
d = 3*exp(-d.^1.3);

figure; imagesc(reshape(d,60,60)');

a = -sum(d);
Ltemp = [ Ltemp, d'; d, a]; %adding the sink node
L{1}.M = -Ltemp;



% %%%%%%%%%%%%%%%%%%  Testing by Denoising Experiment %%%%%%%%%%%%%%%%%%%%%%%
% 
% % cbpdn with graph regularization
mu = .6;
lambda = .16;
lambda_bar = 2;

opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 80;
opt.rho = 100*lambda + 1;
opt.sigma = opt.rho;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;



% Load dictionary
load('CacheData/Dict_12x12.mat');
D = double(D);

opt.Lformat = optl.Lformat;
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

                           
                           


%nonlocal
[Xnl,~]= cbpdnL_lasso_sn(D,sh,L,lambda,lambda_bar,mu,opt);
shrecnl = scnv(D,Xnl);
figure; imagesc(shrecnl); colorbar; title('nonlocal high rec'); colormap(gray);
figure; imagesc(shrecnl+sl); colorbar; title('nonlocal rec'); colormap(gray);
figure; imagesc(sum3(abs(Xnl))); colorbar; title('nonlocal coef');colormap(gray);

[Xnl,~]= cbpdn(D,sh,lambda,opt);
shrecnl = scnv(D,Xnl);
figure; imagesc(shrecnl); colorbar; title(' high rec'); colormap(gray);
figure; imagesc(shrecnl+sl); colorbar; title(' rec'); colormap(gray);
figure; imagesc(sum3(abs(Xnl))); colorbar; title(' coef');colormap(gray);

