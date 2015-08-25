%test the difference between Lpsc and x^2 regularizer; 

load('Euclideaneig1.mat');
imsz = [256,256];
L{1}.ind1 = [1,1];
L{1}.ind2 = imsz;
L{1}.phi = phi;
E = (E./max(E)).^.5;
L{1}.E = E;

s = double(stdimage('lena.grey')) / 255;
s = imresize(s,.5);
[sl,sh] = lowpass(s,4,15);

load('CacheData/Dict_12x12.mat');
D = double(D);


% % cbpdn with graph regularization
mu = .1;
lambda = .15;

opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 100;
opt.rho = 100*lambda + 1;
opt.sigma = opt.rho;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;
opt.Lformat = 'Eig';


scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

%nonlocal
[Xnl,~]= cbpdnL_split(D,sh,L,lambda,mu,opt);
shrecnl = scnv(D,Xnl);
figure; imagesc(shrecnl); colorbar; title('nonlocal high rec'); colormap(gray);
figure; imagesc(shrecnl+sl); colorbar; title('nonlocal rec'); colormap(gray);
figure; imagesc(sum3(abs(Xnl))); colorbar; title('nonlocal coef');colormap(gray);


% x2
opt.HighMemSolve = 1;
[Xx2,~] = cbpdnx2(D,sh,lambda,mu,opt);
shrec = scnv(D,Xx2);
figure; imagesc(shrec); colorbar; title('conv x2 high rec'); colormap(gray);
figure; imagesc(shrec+sl); colorbar; title('conv x2 rec'); colormap(gray);
figure; imagesc(sum3(abs(Xx2))); colorbar; title('conv x2 coef');colormap(gray);

norm(vec(Xnl-Xx2))/norm(vec(Xx2))















