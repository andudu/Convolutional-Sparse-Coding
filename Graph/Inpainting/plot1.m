%plot some stats

load('LenaIp2.mat');
psnr5_nl = squeeze(psnr_rec_nl(1,:,:,:)); 
psnr5_cn = squeeze(psnr_rec_cn(1,:,:,:)); 

l = 2; 

a = []; 
b = [];

for k = 1:length(lambda_all)
    a(k) = max(psnr5_nl(:,k,l));
    b(k) = max(psnr5_cn(:,k,l));    
end

plot(lambda_all, a, 'r', lambda_all, b, 'b'); 
xlabel('lambda'); 
ylabel('psnr') ; 
legend('Laplacian', 'Regular') ; 
title('psnr vs lambda, ndict = 25, noise = .5');

