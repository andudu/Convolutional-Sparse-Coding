%load precomputed graph info

load('data_graph.mat');

% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
sh = sh+.1*randn(size(sh));


%set up parameters
lambda = 0.15;
mu = 100;
opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 500;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;

%optimize
[X,~]= cbpdn_L(D,sh,L,lambda,mu,opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstructing conv
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');
shrec = scnv(D,X);

figure; 
imagesc(shrec);
colorbar;

figure; 
imagesc((sum3(abs(X))));
colorbar;

