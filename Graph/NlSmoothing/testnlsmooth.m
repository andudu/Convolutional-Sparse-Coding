% scrit for testing nlsmooth

%%%%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Coef2.mat');

% full matrix
load([sporco_path,'/Graph/CacheData/Standard/Mat/CosineM1.mat']);
ltemp = L;
L = {};
L{1,1}.M= ltemp{1}.M;
L{1,1}.ind1 = [1,1];
L{1,1}.ind2 = [256,256];


%%%%%%%%%%%%%%%%%%%%%%%%%% Set Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.L1 = 0;





%%%%%%%%%%%%%%%%%%%%%%%%%% Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure; 
imagesc(X);
title('Original');

% figure;
% spy(X);
% title('Original');

mu_all = [.3,.6,1,2,3,4,5,10];
for i = 1:length(mu_all)
    mu = mu_all(i);
    y = nlsmooth(X,mu,L,opt);
    figure;
    imagesc(y);
    title(['Diffused with mu = ',num2str(mu)]);
end



% figure;
% y(y<3e-5) = 0;
% spy(y);
% title(['Diffused with mu = ',num2str(mu)]);