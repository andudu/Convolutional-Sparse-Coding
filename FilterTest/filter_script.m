% gr = [1,-1];
% gr = padarray(gr,[100,100],'post');
% test = dct2(gr);
% test = test.*conj(test);
% figure;
% imagesc(test);
% 
% 
% 
% t = 2*pi*1i*(0:0.01:1);
% figure;
% plot(abs(1-exp(-t)));
% 



% Load test image
s = single(rgbtogrey(stdimage('barbara')))/255;
if isempty(s),
  error('Data required for demo scripts has not been installed.');
end
s = imresize(s,[512,512]);

figure;
imagesc(s);
colormap(gray);

% foo = fftshift(fft2(s));

foo = fft2(s);


a = [size(foo,1)/3,size(foo,2)/3];
foo = foo(1:a(1),1:a(2));
foo = padarray(foo,[256*2,256*2],'post');
c = ifft2(foo);
figure;
imagesc(abs(c));
colormap(gray);

