% generate
global sporco_path;

% load the image
S0 = double(stdimage('lena.grey')) / 255;
S0 = imresize(S0,.5);
sigma = .04;
S0 = S0+sigma*randn(size(S0));


% load clean dictionaries
load('CacheData/Dict_12x12x35.mat');
D0 = double(D);


% load the graph laplacian
load([sporco_path,'/Graph/CacheData/Standard/Mat/CosineM1.mat']);


%%%%%%%%%%%%%%%%%%%%%%%%%%  First Batch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_all = [.1,.2,.25,.3,.35,.4];
mu_all = [.2,.4,.6,.8];
maxit = 300;

disp(['batch 1']);



for i = 1:length(mu_all)
    for j = 1:length(lambda_all)
        disp([num2str(i), ',', num2str(j)]);
        D_init = D0(:,:,1:15);
        lambda = lambda_all(j);
        mu = mu_all(i);
        [D1,X1, D2, X2 ] = exp_dict_4(D_init, S0, L, lambda, mu, maxit);
        fname = ([sporco_path,'/Graph/CacheData/DictNoisyLenaCbpdn1/DictComp',...
            num2str(i),num2str(j),'.mat']);
        save(fname,'D2','X2','lambda','D1','X1','mu','sigma');
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%  Second Batch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_dict = 12:4:35;
disp(['batch 2']);

for i = 1:length(num_dict)
    for j = 1:length(lambda_all)
        disp([num2str(i), ',', num2str(j)]);
        n = num_dict(i);
        D_init = D0(:,:,1:n);
        lambda = lambda_all(j);
        mu = lambda*1.5;
        [D1,X1, D2, X2 ] = exp_dict_4(D_init, S0, L, lambda, mu, maxit);
        fname = ([sporco_path,'/Graph/CacheData/DictNoisyLenaCbpdn2/DictComp',...
            num2str(i),num2str(j),'.mat']);
        save(fname,'D2','X2','lambda','D1','X1','mu','sigma','n');
    end
end

