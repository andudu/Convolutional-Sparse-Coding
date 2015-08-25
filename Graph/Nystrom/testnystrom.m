%script for testing Nystrom
%
% %snoise = randn(256,256)*.1;
% sref = double(stdimage('lena.grey')) / 255;
% sref = imresize(sref,.5);
% sref = sref(60:1:200, 60:1:200);
% %s = sref + .1*randn(size(sref));
% imsz = size(sref);
% psz = [8,8];
% stpsz = [1,1];
% s = padarray(sref,psz-[1,1],'symmetric','post'); %padd only psz-1 elements!
% %[sl,sh] = lowpass(s,1,15);
% % s = sh;


% lena patch
s = double(stdimage('lena.grey')) / 255;
%s = imresize(s,.5);
% sref = s(50:1:60+50-1,160:1:60+160-1);
% s = s + .03*randn(size(s));
[sl,sh] = lowpass(s,7,15);
% s = sh(50:1:60+50-1,160:1:60+160-1);


imsz = size(s);
psz = [12,12];
stpsz = [1,1];
sh = padarray(sh,psz-[1,1],'symmetric','post');

%generate graph
opt.tau = .8;
opt.Laplacian = 'n';
opt.numsample = 500;
opt.Metric = 'Euclidean';
opt.neig = 15;


data = imageblocks(sh,psz,stpsz);
data = reshape(data,size(data,1)*size(data,2),size(data,3));
data = data';
[phi,E] = nystrom(data,opt);
disp('Normalized graph generated');
X1 = [];
for i = 1:opt.neig,
    v = phi(:,i);
    v = reshape(v,imsz)';
    X1(:,:,i) = v;
end
square_plot(X1,{});
figure;
plot(E);

% [phi,E] = nystrom_u(data,opt);
% disp('Unnormalized graph generated');
% X = [];
% for i = 2:opt.neig,
%     v = phi(:,i);
%     v = reshape(v,coefsz)';
%     X2(:,:,i-1) = v;
% end
% square_plot(X2,{});

