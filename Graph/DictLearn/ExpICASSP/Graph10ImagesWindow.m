% script for generating the Graph Laplacian of 10 images

% load the image
ind = [1,5,16,20,23,25,38,46,52];
load('FlickrCC_512_512.mat');
S0 = []; 
S0 = S(:,:,ind); 
S0 = double(S0)/255;
temp = [];
for i = 1:size(S0,3)
    temp(:,:,i) = imresize(S0(:,:,i),.5);
    temp(:,:,i) = temp(:,:,i);% add a little noise
end
S0 = temp;
    
clear S; 

% Filter input images and compute highpass images
npd = 16;
fltlmbd = 4;
[Sl, Sh_all] = lowpass(S0, fltlmbd, npd);

for i = 1:size(S0,3)
    
    disp(num2str(i));
    
    tauc = 1;
    
    tag1 = 'Flicker9Image_WindowKnearest';
    if ~exist([sporco_path,'/Graph/CacheData/',tag1],'dir')
        mkdir([sporco_path,'/Graph/CacheData/',tag1]);
    end
    Sh = Sh_all(:,:,i);
    
    imsz = size(Sh);
    psz = [12,12];
    stpsz = [1,1];
    %generate graph
    optl = {};
    optl.wsz = imsz;
    optl.psz = [12,12];
    optl.neig = [];
    optl.Lformat = 'Sparse';
    optl.Graph.Laplacian = 'n';
    optl.Graph.tau = tauc;
    optl.Graph.Metric = 'Cosine';
    optl.Graph.GraphType = 'KNearest';
    optl.SaveMem = 1;
    optl.Graph.nsz = 14;
    optl.Graph.k = 30;
    
    [L,~] = laplacian_from_image(Sh,optl); 
    
    if ~exist([sporco_path,'/Graph/CacheData/',tag1,'/Mat'],'dir')
        mkdir([sporco_path,'/Graph/CacheData/',tag1,'/Mat']);
    end
    save([sporco_path,'/Graph/CacheData/',tag1,'/Mat', '/',optl.Graph.Metric,'M',num2str(i),'.mat'],'L');       

    disp(['Graph Generated']);     
    
    
    
end