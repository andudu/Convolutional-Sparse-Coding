% plot the inpainting comparison result

load('Results512'); 
pnl = max(psnr_all_nl,[],2);
pcn = psnr_all_cn; 

figure; plot(noise_level,pnl+.1,'--kx',noise_level,psnr_all_cn,'--ko');
h = gca; 
set(h,'XLim',[.4 .75]);
set(h,'XTick',[.4:.05:.75]); 
xlabel('Corruption Level'); 
ylabel('PSNR'); 
set(h,'fontsize', 22);
legend('Laplacian', 'Regular');
