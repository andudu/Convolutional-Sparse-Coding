% script for testing if noise in dictionaries when training with low lambda
% is caused by joint learning too many patches


%load data
load('/nh/u/mathluo/Convolutional-Sparse-Coding/Graph/CacheData/DictLenaCbpdn2/DictComp81.mat');
S = single(stdimage('lena.grey')) / 255;
[~,Sh] = lowpass(S,4,16);


X = X2(:,:,[1:1:26,28:1:34],:);


opt.MaxMainIter = 200;
opt.Verbose = 1;


% with threshold
X(abs(X)<0.06) = 0;

dsz = [size(D2,1),size(D2,2)];
[D,~] = ccmod(X,Sh,dsz,opt);

o.grey = 1;

square_plot(D,o);

save('81ccmod_thre2.mat','D');