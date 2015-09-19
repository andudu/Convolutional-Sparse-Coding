% noise free all out carpet comparison

global sporco_path;

% load the image

imnum = 4; 
ind = [1,5,16,20,23,25,38,46,52];
ind = ind(1:imnum); 
load('FlickrCC_512_512.mat');
S0 = []; 
S0 = S(:,:,ind); 
S0 = single(S0)/255;
temp = [];
sigma = 0.06;
for i = 1:size(S0,3)
    temp(:,:,i) = imresize(S0(:,:,i),.5);
    temp(:,:,i) = temp(:,:,i) + sigma*randn(size(temp(:,:,i)));% add a little noise
end
S0 = temp;
    
clear S; 


% Filter input images and compute highpass images
npd = 16;
fltlmbd = 4;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);

% Load the Laplacians
temp = {};
for i = 1:imnum
   load([sporco_path,'/Graph/CacheData/Flicker9Image_Knearest/Mat/CosineM',num2str(i),'.mat']);
   temp{i,1} = L{1,1};
end
L = temp;


% Dict 3 The real dictionary
% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D = double(D);
D0 = D(:,:,1:1:25); %top 25 dictionaries



%%%%%%%%%%%%%%%%%%%%%%%%%%  First Batch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_all = [.2,.25,.3,.35];
mu_all = [.1,.2,.3,.4];
maxit = 350;



disp(['batch 1 4 Images']);

tag1 = 'DictRegularFlicker_4Image'; % make a directory
if ~exist([sporco_path,'/Graph/CacheData/',tag1],'dir')
    mkdir([sporco_path,'/Graph/CacheData/',tag1]);
end

for i = 1:length(lambda_all)
    
    lambda = lambda_all(i);
    opt = [];
    opt.Verbose = 0;
    opt.MaxMainIter = maxit;
    opt.rho = 50*lambda + 0.5;
    opt.sigma = size(Sh,3);
    opt.AutoRho = 1;
    opt.AutoRhoPeriod = 10;
    opt.AutoSigma = 1;
    opt.AutoSigmaPeriod = 10;
    opt.XRelaxParam = 1.5;
    opt.DRelaxParam = 1.5;
    [D1, X1, ~] = cbpdndliu(D0, Sh, lambda, opt); %regular dictionary
    save([sporco_path,'/Graph/CacheData/',tag1, '/','RegularDict',num2str(i),'.mat'],'D1','X1','lambda');        
    
    for j = 1:length(mu_all)
        opt = [];
        opt.Verbose = 0;
        opt.MaxMainIter = maxit;
        opt.rho = 50*lambda + 0.5;
        opt.sigma = size(Sh,3);
        opt.AutoRho = 1;
        opt.AutoRhoPeriod = 10;
        opt.AutoSigma = 1;
        opt.AutoSigmaPeriod = 10;
        opt.AutoRhoBar = 1;
        opt.AutoRhoBarPeriod = 10;
        opt.Lformat = 'Sparse';
        opt.XRelaxParam = 1.5;
        opt.DRelaxParam = 1.5;
        disp([num2str(i), ',', num2str(j)]);
        mu = mu_all(j);
        [D2, X2,~ ] = cbpdnLdliu(D0, Sh,L,mu,lambda, opt);
        save([sporco_path,'/Graph/CacheData/',tag1, '/','LaplacianDict',num2str(i),num2str(j),'.mat'],'D2','X2','lambda','mu');  
    end
end









