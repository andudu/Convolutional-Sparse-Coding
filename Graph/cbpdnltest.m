% Testing different cbpdnl solver

%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting options for Graph Laplacian
optl = {};
optl.wsz = [60,60];
optl.psz = [8,8];
optl.Lformat = 'Sparse';
optl.neig = 30;
optl.Laplacian = 'n';
optl.Graph.tau = 2;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'Window';
optl.Graph.nsz = [7,7];
optl.Graph.k = [];


imsz = optl.wsz+optl.psz-[1,1] ;
psz = optl.psz;
stpsz = [1,1];


%%%%%%%%%%%%%%%%%%%%%%%%% Generating Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = double(stdimage('lena.grey')) / 255;
s = imresize(s,.5);
s_ref = s(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
sn = s(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
[sl,sh] = lowpass(s,7,15);
s = sh(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
sl = sl(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);


%%%%%%%%%%%%%%%%%%%%%%%%% Generating Laplacian %%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L,scrop] = laplacian_from_image(s,optl);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Testing Cbpdnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % cbpdn with graph regularization
mu = .8;
lambda = .23;

opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 100;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-4;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;
opt.AutoRho = 1;
opt.AutoSigma = 1;



% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('8x8x32');
D = double(D);
opt.Lformat =  optl.Lformat;


tic;
[Xnl1,~] = cbpdnL_split(D,s,L,lambda,mu,opt);
toc;
tic;
[Xnl2,~] = cbpdnL_lasso(D,s,L,lambda,mu,opt);
toc;



disp(['rel 1 norm error',num2str(sum(abs(vec(Xnl1 - Xnl2))) / abs(sum(vec(Xnl2))))]);
disp(['rel max norm error',num2str(max(abs(vec(Xnl1-Xnl2)))/max(vec(Xnl2)))]);







