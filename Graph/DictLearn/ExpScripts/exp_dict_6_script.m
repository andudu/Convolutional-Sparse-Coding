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

lambda_all = [.1,.2,.3,.35,.4];
mu_all = [.1,.2,.3,.4,.5];
num_dict = [6,10,12,15,18,20];
maxit = 300;

disp(['batch 1']);

for i = 1:length(mu_all)
    for j = 1:length(lambda_all)
        for k = 1:length(num_dict)
            disp([num2str(i), ',', num2str(j),',',num2str(k)]);
            numdict = num_dict(k);
            D_init = D0(:,:,1:num_dict(k));
            lambda = lambda_all(j);
            mu = mu_all(i);
            [D1,X1, D2, X2 ] = exp_dict_4(D_init, S0, L, lambda, mu, maxit);
            fname = ([sporco_path,'/Graph/CacheData/DictLenaCbpdnAll_Window/DictComp',...
                num2str(i),num2str(j),num2str(k),'.mat']);
            save(fname,'D2','X2','lambda','D1','X1','mu','numdict');
        end
    end
end




