%script for testing Nystrom

%snoise = randn(256,256)*.1;
sref = double(stdimage('lena.grey')) / 255;
s = sref; %+snoise;
imsz = size(sref);
psz = [12,12];
stpsz = [1,1];
s = padarray(sref,psz,'symmetric','post');
[sl,sh] = lowpass(s,5,15);


%generate graph
opt.tau = 1;
opt.numsample = 500;
opt.Metric = 'Cosine';
opt.neig = 10;


data = imageblocks(sh,psz,stpsz);
data = reshape(data,size(data,1)*size(data,2),size(data,3));
data = data';
[phi,E] = nystrom_n(data,opt);
disp('Normalized graph generated');
X = [];
for i = 2:opt.neig,
    v = phi(:,i);
    v = reshape(v,imsz)';
    X1(:,:,i-1) = v;
end
square_plot(X1,{});

% [phi,E] = nystrom_u(data,opt);
% disp('Unnormalized graph generated');
% X = [];
% for i = 2:opt.neig,
%     v = phi(:,i);
%     v = reshape(v,coefsz)';
%     X2(:,:,i-1) = v;
% end
% square_plot(X2,{});

