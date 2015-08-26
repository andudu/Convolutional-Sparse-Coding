
%generate
global sporco_path;
D0 = zeros(12,12,35, 'single');
D0(4:9,4:9,:) = single(randn(6,6,35));


load('CacheData/Dict_12x12x35.mat');
D0_oe = D;

%%%%%%%%%%%%%%%%%%%%%%%%%%  First Batch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_all = [.02,.05,.1,.2,.3];
num_dict = 6:4:35;

disp(['batch 1']);

for i = 1:length(num_dict)
    for j = 1:length(lambda_all)
        disp([num2str(i), ',', num2str(j)]);
        n = num_dict(i);
        D_init = D0(:,:,1:n);
        lambda = lambda_all(j);
        imind = 1;
        maxit = 300;
        [D2,X2,Aind2] = exp_dict_3(D_init,lambda,imind,maxit);
        fname = ([sporco_path,'/Graph/CacheData/DictLenaCbpdn1/DictComp',...
            num2str(i),num2str(j),'.mat']);
        save(fname,'D2','Aind2','X2','lambda','imind','n');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%  Second Batch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lambda_all = [.02,.05,.1,.2,.3];
num_dict = 6:4:35;

disp(['batch 2']);

for i = 1:length(num_dict)
    for j = 1:length(lambda_all)
        disp([num2str(i), ',', num2str(j)]);
        n = num_dict(i);
        D_init = D0_oe(:,:,1:n);
        lambda = lambda_all(j);
        imind = 1;
        maxit = 300;
        [D2,X2,Aind2] = exp_dict_3(D_init,lambda,imind,maxit);
        fname = ([sporco_path,'/Graph/CacheData/DictLenaCbpdn2/DictComp',...
            num2str(i),num2str(j),'.mat']);
        save(fname,'D2','Aind2','X2','lambda','imind','n');
    end
end

