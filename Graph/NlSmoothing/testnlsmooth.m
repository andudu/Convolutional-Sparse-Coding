% scrit for testing nlsmooth

%%%%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Coef2.mat');

% full matrix
load(['/nh/u/mathluo/Convolutional-Sparse-Coding/Graph/CacheData/Lena_WindowKnearest_varn/Mat/CosineM1.mat']);
ltemp = L;
L1 = {};
L1{1,1}.M= ltemp{1}.M;
L1{1,1}.ind1 = [1,1];
L1{1,1}.ind2 = [256,256];


load(['/nh/u/mathluo/Convolutional-Sparse-Coding/Graph/CacheData/Lena_WindowKnearest_varn/Mat/CosineM4.mat']);
ltemp = L;
L2 = {};
L2{1,1}.M= ltemp{1}.M;
L2{1,1}.ind1 = [1,1];
L2{1,1}.ind2 = [256,256];




%%%%%%%%%%%%%%%%%%%%%%%%%% Set Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.L1 = 0;





%%%%%%%%%%%%%%%%%%%%%%%%%% Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure; 
imagesc(X);
title('Original');

mu_all = [.3,.6,1,2,3,4,5,10];

for i = 1:length(mu_all)
    figure;
    mu = mu_all(i);
    y = nlsmooth(X,mu,L1,opt);    
    subplot(1,2,1);
    imagesc(y);
    title(['Diffuse k small with mu = ',num2str(mu)]);
    
    y = nlsmooth(X,mu,L2,opt);    
    subplot(1,2,2);
    imagesc(y);
    title(['Diffuse k big with mu = ',num2str(mu)]);    
   
end

