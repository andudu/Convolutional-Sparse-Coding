%script for testing Nystrom

snoise = randn(256,256)*.1;
sref = double(stdimage('lena.grey')) / 255;
sref = imresize(sref,.5);
s = sref+snoise;
[sl,sh] = lowpass(s,5,15);


%generate graph
opt.tau = 1;
opt.numsample = 200;
opt.Metric = 'Cosine';
opt.neig = 10;
imsz = [256,256];
psz = [12,12];
stpsz = [1,1];
coefsz = imsz - psz + [1,1];
data = imageblocks(sh,psz,stpsz);
data = reshape(data,size(data,1)*size(data,2),size(data,3));
data = data';
[phi,E] = nystrom_n(data,opt);
disp('Normalized graph generated');
X = [];
for i = 2:opt.neig,
    v = phi(:,i);
    v = reshape(v,coefsz)';
    X1(:,:,i-1) = v;
end
square_plot(X1,{});

[phi,E] = nystrom_u(data,opt);
disp('Unnormalized graph generated');
X = [];
for i = 2:opt.neig,
    v = phi(:,i);
    v = reshape(v,coefsz)';
    X2(:,:,i-1) = v;
end
square_plot(X2,{});

