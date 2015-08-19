% script for generating L1 Weights

% only perfect case right now. 

% testing if reweighing L1 is feasible for removing white noise


%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load  Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load a saved noise
load('CacheData/stdnoise.mat');
sref = double(stdimage('lena.grey')) / 255;
sref = imresize(sref,.5);
[slref,shref] = lowpass(sref,5,15);
s = sref+r_noise;
[sl,sh] = lowpass(s,5,15);


% Load nice dictionary
load('CacheData/Dict_12x12.mat');




%%%%%%%%%%%%%%%%%%%%%%% Perfect Weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mu = .5;
% lambda = 0.2;
% opt = {};
% opt.Verbose = 1;
% opt.MaxMainIter = 50;
% opt.rho = 10;
% opt.RelStopTol = 1e-3;
% opt.AuxVarObj = 0;
% opt.HighMemSolve = 1;
% opt.L1Weight = 1;
% 
% [X,~]= cbpdngr(D,shref,lambda,mu,opt);
% 
% square_plot(X,{});
% 
% Weight = [];
% 
% for i = 1:size(X,3)
%     temp = abs(X(:,:,i));
%     d = max(vec(temp));
%     temp = temp/d;
%     Weight(:,:,i) = exp(-1.8*temp);    
% end
% 
% square_plot(Weight,{});


%%%%%%%%%%%%%%%%%%%%%%%%% Graph Weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
























































