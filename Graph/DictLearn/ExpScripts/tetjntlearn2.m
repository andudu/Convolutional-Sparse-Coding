% script for testing if noise in dictionaries when training with low lambda
% is caused by joint learning too many patches


%load data

load([sporco_path,'/Graph/CacheData/DictLenaCbpdn2/DictComp53.mat']);
S = single(stdimage('lena.grey')) / 255;
[~,Sh] = lowpass(S,4,16);

X = X2;
%X = X2(:,:,[1:1:26,28:1:34],:); %getting rid of dirac delta

%add some noise to X
ind = randi(256, 2,328); 
for k = 1:size(ind,2)
    X(ind(1,k),ind(2,k),:)  = X(ind(1,k),ind(2,k),:) + randn(1,1,22)*.1;
end

opt.MaxMainIter = 200;
opt.Verbose = 1;
opt.AutoSigma = 1;

% with threshold
%X(abs(X)<0.06) = 0;

dsz = [size(D2,1),size(D2,2)];
[D,~] = ccmod(X,Sh,dsz,opt);

o.grey = 1;

square_plot(D,o);

save('53pertccmod.mat','D');