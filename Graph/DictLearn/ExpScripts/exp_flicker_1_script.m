% flicker 5 image test with k-nearest neighbor graph

global sporco_path;

% load the image
ind = [  28 ,   34 ,   42 ,   61 ,   99 ];
load('Flicker1_512_split.mat');
S0 = S(:,:,ind); 
clear S; 


% load clean dictionaries
load('CacheData/Dict_12x12x35.mat');
D0 = double(D);


% load the graph laplacian
Ltemp = {}; 
for i = 1:length(ind)
    load([sporco_path,'/Graph/CacheData/Flicker_KNearest_Cosine/Mat/CosineM',num2str(ind(i)),'.mat']);
    Ltemp{i,1}.M = L{1,1}.M; 
    Ltemp{i,1}.ind1 = L{1,1}.ind1; 
    Ltemp{i,1}.ind2 = L{1,1}.ind2;  
end
L = Ltemp; 
clear Ltemp;

%%%%%%%%%%%%%%%%%%%%%%%%%%  First Batch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_all = .1:.05:.3;
mu_all = .15:.05:.25;
num_dict = 15:5:35;
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
            fname = ([sporco_path,'/Graph/CacheData/DictFlicker5All/DictComp',...
                num2str(i),num2str(j),num2str(k),'.mat']);
            save(fname,'D2','X2','lambda','D1','X1','mu','numdict','ind');
        end
    end
end


