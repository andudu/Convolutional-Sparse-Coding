%plot some stats

load('LenaIp2.mat');
psnr5_nl = squeeze(psnr_rec_nl(1,:,:,:)); 
psnr5_cn = squeeze(psnr_rec_cn(1,:,:,:)); 

 

a = []; 
b = [];

for k = 1:length(numdict)
    a(k) = max(vec(psnr5_nl(:,:,k)));
    b(k) = max(vec(psnr5_cn(:,:,k)));    
end

plot(numdict, a, 'r', numdict, b, 'b'); 
xlabel('lambda'); 
ylabel('psnr') ; 
legend('Laplacian', 'Regular') ; 
title('psnr vs lambda, ndict = 25, noise = .5');

