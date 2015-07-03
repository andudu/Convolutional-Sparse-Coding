% Script demonstrating usage of the cbpdndliu function.
%
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2015-04-09
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'Copyright' and 'License' files
% distributed with the library.


% Training images
S0 = zeros(512, 512, 5, 'single');
S0(:,:,1) = single(stdimage('lena.grey')) / 255;
S0(:,:,2) = single(stdimage('barbara.grey')) / 255;
S0(:,:,3) = single(stdimage('kiel.grey')) / 255;
S0(:,:,4) = single(rgb2gray(stdimage('mandrill'))) / 255;
tmp = single(stdimage('man.grey')) / 255;
S0(:,:,5) = tmp(101:612, 101:612);


%Reduce images size to speed up demo script
tmp = zeros(256, 256, 5, 'single');
for k = 1:size(S0,3),
  tmp(:,:,k) = imresize(S0(:,:,k), 0.5);
end
S0 = tmp;


% Filter input images and compute highpass images
npd = 16;
fltlmbd = 3.5;
%[Sl, Sh] = lowpass(S0, fltlmbd, npd);


% Construct initial dictionary
num_dict_high = 12; %smaller dictionary
num_dict_mid = 1;
num_dict_low = 1;
num_dict = num_dict_high+num_dict_low+num_dict_mid;
small_filter_sz = 8;
mid_filter_sz = 16;
big_filter_sz = 32;


D0 = zeros(big_filter_sz,big_filter_sz,num_dict, 'single');

%initialization for high,mid,low
D0(small_filter_sz/4+1:small_filter_sz/4*3,small_filter_sz/4+1:small_filter_sz/4*3 ...
,1:num_dict_high) = single(randn(small_filter_sz/2,small_filter_sz/2,num_dict_high));

D0(mid_filter_sz/4+1:mid_filter_sz/4*3,mid_filter_sz/4+1:mid_filter_sz/4*3 ...
,num_dict_high+1:num_dict_high+num_dict_mid) = single(randn(mid_filter_sz/2,mid_filter_sz/2,num_dict_mid));

D0(big_filter_sz/4+1:big_filter_sz/4*3,big_filter_sz/4+1:big_filter_sz/4*3 ...
,num_dict_high+num_dict_mid+1:end) = single(randn(big_filter_sz/2,big_filter_sz/2,num_dict_low));

% Set up cbpdndliu parameters
lambda = 0.2;
mu1= 100; %grd on Dict
mu2 = 0;    %grd on Coef
wl1 = [1.3,2.2];  
wgrd = [.45,1];

opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(S0,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 6;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 6;
opt.XRelaxParam = 1.3;
opt.DRelaxParam = 1.3;
opt.GrCoef = 1;

opt.L1Weight = reshape([ones(1,num_dict_high),wl1],1,1,num_dict);
opt.DgrdWeight = reshape([zeros(1,num_dict_high),wgrd],1,1,num_dict);
%opt.XgrdWeight = opt.DgrdWeight;
opt.DictFilterSizes = [repmat([small_filter_sz;small_filter_sz],1,num_dict_high),...
    repmat([mid_filter_sz;mid_filter_sz],1,num_dict_mid),...
    repmat([big_filter_sz;big_filter_sz],1,num_dict_low)];

% Do dictionary learning
[D, X, optinf] = cbpdndliu_grd(D0, S0, mu1,mu2, lambda, opt);


% Plot functional value evolution
figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');

% Display learned dictionary
opt1.grey = 1;
opt1.unifscale = 0;
h = cell(1,4);
h{1} = square_plot(D,opt1);

%plot the "high freq" plus "low freq" for all stuff
DX = zeros(size(X,1),size(X,2),size(X,4));
for k = 1:5
    DX(:,:,k) = convsum(D,X(:,:,:,k),1:1:num_dict_high);
end
opt1.grey = 1;
opt1.unifscale = 0;
h{2} = square_plot(DX,opt1);

%mid
DX = zeros(size(X,1),size(X,2),size(X,4));
for k = 1:5
    DX(:,:,k) = convsum(D,X(:,:,:,k),num_dict_high+1:num_dict_high+num_dict_mid);
end
opt1.grey = 1;
opt1.unifscale = 0;
h{3} = square_plot(DX,opt1);



%low
DX = zeros(size(X,1),size(X,2),size(X,4));
for k = 1:5
    DX(:,:,k) = convsum(D,X(:,:,:,k),num_dict_high+num_dict_mid+1:num_dict);
end
opt1.grey = 1;
opt1.unifscale = 0;
h{4} = square_plot(DX,opt1);


%saving the work
folder_tag = 'ThreeLayerDict256';
tag = 'D';
if ~exist(strcat(sporco_path,'/Results/DictLearn/',folder_tag),'dir')
    mkdir(strcat(sporco_path,'/Results/DictLearn/',folder_tag));
end
name = {'Dict','LowRec','MidRec','HighRec'};
for i = 1:4
    saveas(h{i},strcat(sporco_path,'/Results/DictLearn/',folder_tag,'/',...
        tag,name{i}),'fig');
end

par.mu1 = mu1;
par.mu2 = mu2;
par.lambda = lambda;
parl.wl1 = wl1;
par.wgrd = wgrd;
par.fsize = opt.DictFilterSizes;
save(strcat(sporco_path,'/Results/DictLearn/',folder_tag,'/Dict',tag,'.mat'),'D','par');













