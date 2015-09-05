% Graph Generation Script 5
% Script for generating a series of graph from the Lena image
% unnormalized Graph varying k


%%%%%%%%%%%%%%%%%%%%%%%%%%% Batch 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('Flicker1_512_split.mat');

tag1 = 'Flicker_KNearest_Cosine';
if ~exist([sporco_path,'/Graph/CacheData/',tag1],'dir')
    mkdir([sporco_path,'/Graph/CacheData/',tag1]);
end


lambda = 4;
tauc = 1;


% Building the Graph

for i = 1:100
    i
    S0 = S(:,:,i);
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
    optl.Graph.k = 30;
    
    disp('generating Cosine graph');
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
