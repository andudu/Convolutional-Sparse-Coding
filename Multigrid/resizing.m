% Load test image
s = single(rgbtogrey(stdimage('barbara')))/255;
if isempty(s),
  error('Data required for demo scripts has not been installed.');
end
s1 = imresize(s,[512,512]);

c = 4;
s2 = imresize(s1,1/c);
s3 = imresize(s2,c);

figure;
imagesc(s2);
figure;
imagesc(s3);


