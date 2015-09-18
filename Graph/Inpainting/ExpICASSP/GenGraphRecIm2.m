% Generate some graphs from some reconstructed images

% load the standard image
s_ref = imresize(double(stdimage('lena.grey'))/255,.5);

% load the dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D = double(D);
temp = D;
D(:,:,1) = zeros(12,12); 
D(1,1,1) = 1; 
D(:,:,2:end+1) = temp; 


% set parameters
noise_level = [.4,.5,.6,.65,.7,.75]; 
mu1_best = [.5,.5,.5,.5,.5,.5];  % mu1 for the lower frequency
lambda_best = [.005,.01,.01,.01,.02,.02]; 
maxiter = 500;



for i = 1:length(noise_level) %different noise level
    disp([num2str(i)]);
    mu1 = mu1_best(i);
    lambda = lambda_best(i);
    
    load(['CorImNoise256',num2str(i),'.mat']);
    
    opt = [];
    opt.L1Weight = ones(size(S_c,1),size(S_c,2),size(D,3));
    opt.L1Weight(:,:,1) = 10*ones(size(S_c,1),size(S_c,2));
    for t = 1:size(ind,2)
        opt.L1Weight(ind(1,t),ind(2,t),1) = 0;
    end
    opt.Verbose = 0;
    opt.MaxMainIter = maxiter;
    opt.RelStopTol = 1e-3;
    opt.rho = 50*lambda + 1;
    opt.AutoRho = 1;
    opt.AutoRhoPeriod = 1;
    opt.RhoRsdlTarget = 0.1;
    opt.RelaxParam = 1.8;
    
    
    [X, Z, ~] = cbpdnlc(D, S_c, lambda, mu1, opt);
    
    
    DX = ifft2(fft2(D, size(X,1), size(X,2)) .* fft2(X), 'symmetric');
    Sh = sum(DX(:,:,2:end), 3);
    s1 = Sh + Z; %final reconstruction
   
    disp(['Reconstruction Complete']);
    disp(['psnr = ', num2str(psnr(s1,s_ref))]);
    disp(['ssim = ', num2str(ssim(s1,s_ref))]);
    
    
    % Graph Generation
    
    tauc = 1;
    
    tag1 = 'Lena256_Knearest_Ip';
    if ~exist([sporco_path,'/Graph/CacheData/',tag1],'dir')
        mkdir([sporco_path,'/Graph/CacheData/',tag1]);
    end
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
    optl.Graph.nsz = [];
    optl.Graph.k = 30;
    
    [L,~] = laplacian_from_image(Sh,optl); 
    
    if ~exist([sporco_path,'/Graph/CacheData/',tag1,'/Mat'],'dir')
        mkdir([sporco_path,'/Graph/CacheData/',tag1,'/Mat']);
    end
    save([sporco_path,'/Graph/CacheData/',tag1,'/Mat', '/',optl.Graph.Metric,'M',num2str(i),'.mat'],'L');       

    disp(['Graph Generated']); 
    
    
end
