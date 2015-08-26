
%generate
global sporco_path;
D0 = zeros(12,12,30, 'single');
D0(4:9,4:9,:) = single(randn(6,6,30));

%%%%%%%%%%%%%%%%%%%%%%%%%%  First Batch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = .2;
mu = .18;
num_dict = [6:3:27];
ind_all = num2cell(1:1:100);

disp(['batch 1']);

for i = 6:length(num_dict)
    for j = 1:length(ind_all)
        disp([num2str(i), ',', num2str(j)]);
        n = num_dict(i);
        D_init = D0(:,:,1:n);
        imind = ind_all{j};
        maxit = 300;
        [D1,D2,Aind1,Aind2] = exp_dict_1(D_init,lambda,mu,imind,maxit);
        fname = ([sporco_path,'/Graph/CacheData/DictFlicker1/Dictcomp1',...
            num2str(i),num2str(j),'.mat']);
        save(fname,'D1','D2','Aind1','Aind2','lambda','mu','imind','n');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%  Second Batch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = .2;
mu = .18;
num_dict = [6:3:27];
ind_all = {};
for i = 1:40
   ind_all{i} = randi(100,1,2); 
end

disp(['batch 2']);

for i = 1:length(num_dict)
    for j = 1:length(ind_all)
        disp([num2str(i), ',', num2str(j)]);
        n = num_dict(i);
        D_init = D0(:,:,1:n);
        imind = ind_all{j};
        maxit = 300;
        [D1,D2,Aind1,Aind2] = exp_dict_1(D_init,lambda,mu,imind,maxit);
        fname = ([sporco_path,'/Graph/CacheData/DictFlicker1/Dictcomp2',...
            num2str(i),num2str(j),'.mat']);
        save(fname,'D1','D2','Aind1','Aind2','lambda','mu','imind','n');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Third Batch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = .2;
mu = .18;
num_dict = [6:2:24];
ind_all = {[1,2,3,4],[1,2,4,5],[1,2,3,5],[2,3,4,5]}; %cross  validate

disp(['batch 3']);

for i = 1:length(num_dict)
    for j = 1:length(ind_all)
        disp([num2str(i), ',', num2str(j)]);
        n = num_dict(i);
        D_init = D0(:,:,1:n);
        imind = ind_all{j};
        maxit = 300;
        [D1,D2,Aind1,Aind2] = exp_dict_1(D_init,lambda,mu,imind,maxit);
        fname = ([sporco_path,'/Graph/CacheData/DictFlicker1/Dictcomp3',...
            num2str(i),num2str(j),'.mat']);
        save(fname,'D1','D2','Aind1','Aind2','lambda','mu','imind','n');
    end
end

