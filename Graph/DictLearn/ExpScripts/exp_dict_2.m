%experiment 2 for dictionary learning
function [D1,D2,Aind1, Aind2] = exp_dict_2(D0,lambda,mu, imind,maxit)

%%%%%%%%%%%%%%%%%%%%%%%%%% Load the Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load([sporco_path,'Graph/CacheData/']);
global sporco_path;
% Flicker Training images
load('Flicker1_512_split.mat');
S0 = S;
S0 = double(S0)/255;

% Reduce images size to speed up demo script
tmp = [];
kk = 1;
for k = imind,
  tmp(:,:,kk) = S0(:,:,k);
  kk = kk+1;
end
S0 = tmp;


ii = 1;
for i = imind
    load([sporco_path,'/Graph/CacheData/Flicker/Eig/Cosineeig',num2str(i),'.mat']);
    L{ii,1}.phi = phi; 
    L{ii,1}.E = (E./max(E)).^0.75;
    L{ii,1}.ind1 = [1,1];
    L{ii,1}.ind2 = [size(S0,1),size(S0,2)];
    ii = ii+1;
end



%%%%%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Filter input images and compute highpass images
npd = 16;
fltlmbd = 4;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);
%Sh = Sh + .05*randn(size(Sh));

% Construct initial dictionary
% numdict = 15;
% D0 = zeros(12,12,numdict, 'single');
% D0(4:9,4:9,:) = single(randn(6,6,numdict));
% numdict = size(D0,3);

% Set up cbpdndliu parameters
%lambda = 0.2;
%mu = lambda/2;
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
opt.AutoRhoBar = 1;
opt.AutoRhoBarPeriod = 10;
opt.Lformat = 'Eig';



%%%%%%%%%%%%%%%%%%%%%%%% Dictionary Learning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do dictionary learning
L = {};
L{1,1}.phi = phi;
L{1,1}.E = E;
L{1,1}.ind1 = [1,1];
L{1,1}.ind2 = [size(S0,1),size(S0,2)];

[D1, X1, ~] = cbpdnLdliu(D0, Sh,L, lambda,mu, opt);
opt = rmfield(opt,{'AutoRhoBar','AutoRhoBarPeriod','Lformat'});
[D2, X2, ~] = cbpdndliu(D0, Sh, lambda, opt);

%%%%%%%%%%%%%%%%%%%%%%%% Recording Unused Dicts %%%%%%%%%%%%%%%%%%%%%%%%%%%
Aind1 = [];
Aind2 = [];
for i = 1: size(D0,3)
    if(nnz(X1(:,:,i,:))>=20)
       Aind1 = [Aind1,i]; 
    end
    if(nnz(X2(:,:,i,:))>=20)
       Aind2 = [Aind2,i]; 
    end      
end


end

