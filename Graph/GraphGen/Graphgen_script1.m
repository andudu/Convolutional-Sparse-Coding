% Graph Generation Script 1
% Script for generating a series of graph from the Lena image
% normalized Graph Laplacian
% varing k in K-nearest neighbor graphs


%%%%%%%%%%%%%%%%%%%%%%%%%%% Batch 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Standard Training images
S0 = double(stdimage('lena.grey')) / 255;
S0 = imresize(S0,.5);

tag1 = 'Lena_Knearest_vark';
if ~exist([sporco_path,'/Graph/CacheData/',tag1],'dir')
    mkdir([sporco_path,'/Graph/CacheData/',tag1]);
end

k_all = [10,20,30,40];
lambda = 4;
tauc = 1;
taue = 5;

% Building the Graph

for i = 1:length(k_all)
    i
    [~,sh] = lowpass(S0,lambda,15);
    imsz = size(sh);
    psz = [12,12];
    stpsz = [1,1];
    %generate graph
    optl = {};
    optl.wsz = imsz;
    optl.psz = [12,12];
    optl.neig = 100;
    optl.Lformat = 'Sparse';
    optl.Graph.Laplacian = 'n';
    optl.Graph.tau = tauc;
    optl.Graph.Metric = 'Cosine';
    optl.Graph.GraphType = 'KNearest';
    optl.SaveMem = 1;
    optl.Graph.nsz = [];
    optl.Graph.k = k_all(i);
    
    disp('generating Cosine graph');
    [L,~] = laplacian_from_image(sh,optl);        
    if ~exist([sporco_path,'/Graph/CacheData/',tag1,'/Mat'],'dir')
        mkdir([sporco_path,'/Graph/CacheData/',tag1,'/Mat']);
    end
    save([sporco_path,'/Graph/CacheData/',tag1,'/Mat', '/',optl.Graph.Metric,'M',num2str(i),'.mat'],'L');    
    
    
    disp('generating Euclidean graph');  
    optl.Graph.tau = taue;
    optl.Graph.Metric = 'Euclidean';
    [L,~] = laplacian_from_image(sh,optl);         
    if ~exist([sporco_path,'/Graph/CacheData/',tag1,'/Mat'],'dir')
        mkdir([sporco_path,'/Graph/CacheData/',tag1,'/Mat']);
    end
    save([sporco_path,'/Graph/CacheData/',tag1,'/Mat', '/',optl.Graph.Metric,'M',num2str(i),'.mat'],'L');       
    
    clear L;
    
end

%save data
imsz = [size(S0,1),size(S0,2)];
save([sporco_path,'/Graph/CacheData/',tag1,'/param.mat'],'lambda','imsz','psz','optl','k_all','tauc','taue');






