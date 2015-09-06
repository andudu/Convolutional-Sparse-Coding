% noise free all out carpet comparison for windowed graph

global sporco_path;

% load the image
S0 = double(stdimage('lena.grey')) / 255;
S0 = imresize(S0,.5);
%sigma = .04;
%S0 = S0+sigma*randn(size(S0));


% load clean dictionaries
load('CacheData/Dict_12x12x35.mat');
D0 = double(D);


% load the graph laplacian
load([sporco_path,'/Graph/CacheData/Lena_WindowKnearest_varn/Mat/CosineM2.mat']);


%%%%%%%%%%%%%%%%%%%%%%%%%%  First Batch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_all = .1:.02:.3;
mu_all = .15:.02:.35;
maxit = 300;

disp(['batch 1']);

for i = 1:length(mu_all)
    for j = 1:length(lambda_all)
            disp([num2str(i), ',', num2str(j)]);
            numdict = 25;
            D_init = D0(:,:,1:numdict);
            lambda = lambda_all(j);
            mu = mu_all(i);
            [D1,X1, D2, X2 ] = exp_dict_4(D_init, S0, L, lambda, mu, maxit);
            fname = ([sporco_path,'/Graph/CacheData/DictLenaCbpdnAll_WindowFine/DictComp',...
                num2str(i),num2str(j),'.mat']);
            save(fname,'D2','X2','lambda','D1','X1','mu','numdict');
    end
end




