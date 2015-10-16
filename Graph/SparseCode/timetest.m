% script for doing a speed test of various algorithms. 
% tests the speed of these four combinations 
% eig-admm 
% eig-lasso
% dir-lasso
% dir-admm
% original cbpdn

% load the image
s = double(stdimage('lena.grey')) / 255;
s = imresize(s,.5); 
[Sl,Sh] = lowpass(s,5,14); 

% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D = double(D); 

eig_time = [];
M_time = []; 
cbpdn_time = []; 

lambda = .04; 
mu = .03; 
opt.Verbose = 0;
opt.MaxMainIter = 4;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.AutoRhoPeriod = 10;



% %lasso-eig
% for i = 1:4
%     %load the eigenvectors
%     opt.Lformat = 'Eig';
%     load(strcat(sporco_path,'/Graph/CacheData/Lena_timetest/Eig',num2str(i),'.mat'),'L','numeig');
%     tic;
%     [Xnl,~] = cbpdnL_lasso(D,Sh,L,lambda,mu,opt);
%     eig_time(i,1) = toc;
% end


% %admm-eig 
% for i = 1:4
%     %load the eigenvectors
%     opt.Lformat = 'Eig'; 
%     load(strcat(sporco_path,'/Graph/CacheData/Lena_timetest/Eig',num2str(i),'.mat'));     
%     tic; 
%     [Xnl,~] = cbpdnL_split(D,Sh,L,lambda,mu,opt);
%     eig_time(i,2) = toc;
% end

%lasso-Full
for i = 1:4
    %load the eigenvectors
    opt.Lformat = 'Mat';
    load(strcat(sporco_path,'/Graph/CacheData/Lena_timetest/EuclideanM',num2str(i),'.mat'));     
    tic; 
    [Xnl,~] = cbpdnL_lasso(D,Sh,L,lambda,mu,opt);
    M_time(i,1) = toc;
end

%admm-Full
for i = 1:4
    %load the eigenvectors
    opt.Lformat = 'Mat';
    load(strcat(sporco_path,'/Graph/CacheData/Lena_timetest/EuclideanM',num2str(i),'.mat'));     
    tic; 
    [Xnl,~] = cbpdnL_split(D,Sh,L,lambda,mu,opt);
    M_time(i,2) = toc;
end



tic; 
[Xnl,~] = cbpdn(D,Sh,lambda,opt); 
cbpdn_time = toc; 

save('time_results.mat','eig_time','M_time','cbpdn_time'); 







  

