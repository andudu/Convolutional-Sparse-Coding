% script for parameter search for regular dictionary learning


% load the standard image
s_ref = single(stdimage('lena.grey'))/255;

% load the dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x64');
temp = D;
D(:,:,1) = zeros(12,12); 
D(1,1,1) = 1; 
D(:,:,2:end+1) = temp; 


% set parameters
noise_level = [.4,.5,.6,.65,.7,.75]; 
mu1_all = [.5,.75,1];  % mu for the lower frequency
lambda_all = [0.006, 0.008, 0.01,0.012,0.014, 0.0016]; 


psnr_all = []; 


for i = 4:length(noise_level) %different noise level
    load(['CorImNoise',num2str(i),'.mat']);
    opt = [];
    
    opt.L1Weight = ones(size(S_c,1),size(S_c,2),size(D,3));
    opt.L1Weight(:,:,1) = 10*ones(size(S_c,1),size(S_c,2));
    
    for t = 1:size(ind,2)
    opt.L1Weight(ind(1,t),ind(2,t),1) = 0;
    end
    
    for j = 1:length(mu1_all)
        mu1 = mu1_all(j);
        for k = 1:length(lambda_all)
            
            disp([num2str(i), num2str(j), num2str(k)]);
            
            lambda = lambda_all(k); 
            opt.Verbose = 0;
            opt.MaxMainIter = 500;
            opt.RelStopTol = 1e-3;
            opt.rho = 50*lambda + 1;
            opt.AutoRho = 1;
            opt.AutoRhoPeriod = 1;
            opt.RhoRsdlTarget = 0.1;
            opt.RelaxParam = 1.8;
           
            
            [X, Z, ~] = cbpdnlc(D, S_c, lambda, mu1, opt);
            
            DX = ifft2(fft2(D, size(X,1), size(X,2)) .* fft2(X), 'symmetric');
            s1 = sum(DX(:,:,2:end), 3) + Z; %final reconstruction
            
            psnr_all(i,j,k) = psnr(s1,s_ref);
            disp(['psnr = ', num2str(psnr_all(i,j,k))]);
            
        end
    end
end

save('Results512_64.mat','psnr_all', 'noise_level','mu1_all','lambda_all'); 
