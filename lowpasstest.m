% %script for choosing a suitable lowpass parameter
% 
% %script 1
% %%%%%%%%%%%%%%%% Data Construction%%%%%%%%%%%%%%%%%
% 
% %construct two images with a line
% s1 = zeros(256,256);
% s1(:, 128:129)= 1;% one pixel wide line
% 
% s2 = zeros(128,128);
% s2(:, 64:65)= 1;% one pixel wide line
% 
% s3 = imresize(s1,.5);
% 
% 
% %%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lambda = 3;
% pad = 15;
% 
% [sl1,sh1] = lowpass(s1,lambda,pad);
% [sl2,sh2] = lowpass(s2,lambda,pad);
% [sl3,sh3] = lowpass(s3,lambda,pad);
% sh2_resz = imresize(sh2,2);
% 
% %%%%%%%%%%%%%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure;
% imagesc(sh1);
% colorbar;
% title('highpass of large image');
% 
% figure;
% imagesc(sh2);
% colorbar;
% title('highpass of small image');
% 
% figure;
% imagesc(sh2_resz);
% colorbar;
% title('upscaled highpass of small image');
% 
% figure;
% imagesc(sh3);
% colorbar;
% title('downscaled highpass');
% 
% 
% % % script 2
% line_width = 10;
% imsz = 256;
% s = zeros(imsz,imsz);
% mid = ceil(imsz/2);
% s(:,mid:mid+line_width)=1;
% 
% for i = 2:2:
%     [~,sh] = lowpass(s,i,15);
%     figure;
%     imagesc(sh);
%     colorbar;
%     title(['lambda = ', num2str(i)]);
% end

% script 3
% img0 = single(stdimage('lena.grey'));
% img0 = img0 / 255;

img0 = randn(256,256)*.1;
sl_mean = meanfilter(img0,1,8);
sh_mean = img0-sl_mean;

figure;
imagesc(sl_mean);
title('mean low freq');
colorbar;
figure;
imagesc(sh_mean);
title('mean high freq');
colorbar;
for i = 20:4:40
    [sl,sh] = lowpass(img0,i,10);
    figure;
    imagesc(sl);
    title(['low lambda = ', num2str(i)]);
    colorbar;
    figure;
    imagesc(sh);
    title(['high lambda = ', num2str(i)]);  
    colorbar;
end

































