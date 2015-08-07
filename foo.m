load('lightning.mat');
a = double(imresize(x,[512,512]));


[al,ah] = lowpass(a,14,13);

figure;
imagesc(al);
colormap(gray);
axis off;
saveas(gcf,'lightninglow','png');

figure;
imagesc(ah);
colormap(gray);
axis off;
saveas(gcf,'lightninghigh','png');