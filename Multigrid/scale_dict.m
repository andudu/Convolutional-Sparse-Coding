%Script for Training Smooth Dictionaries

% % Training images from Repository
crpsz = 512;

%saving the work
folder_tag = 'Test';
tag = 'test';


image_names = { 'airplane.png' ,   'lena.grey.png', 'boats.png' ,'house.png' ,...
    'monarch.png','barbara.grey.png',  'bridge.grey.png',   'kiel.grey.png',...
    'man.grey.png', 'peppers.png', 'goldhill.png',...
    'mandrill.png' };
image_index = 6;
%image_index =  [2,6,7,8,9,10,11];
%image_index = [2,6]; %barbara
image_num = length(image_index);

S0 = zeros(crpsz, crpsz, image_num , 'single');

ind = 1;
for i = image_index
    tmp = stdimage(image_names{i});
    if ndims(tmp) == 3
        tmp = rgb2gray(tmp);
    end
    tmp = imresize(tmp,[crpsz,crpsz]);
    S0(:,:,ind) = single(tmp)/255.0;
    ind = ind+1;
end
clear tmp tmp1; 


% % Std 5 Training images
% S0 = zeros(512, 512, 5, 'single');
% S0(:,:,1) = single(stdimage('lena.grey')) / 255;
% S0(:,:,2) = single(stdimage('barbara.grey')) / 255;
% S0(:,:,3) = single(stdimage('kiel.grey')) / 255;
% S0(:,:,4) = single(rgb2gray(stdimage('mandrill'))) / 255;
% tmp = single(stdimage('man.grey')) / 255;
% S0(:,:,5) = tmp(101:612, 101:612);
% %Reduce images size to speed up demo script
% tmp = zeros(256, 256, 5, 'single');
% for k = 1:size(S0,3),
%   tmp(:,:,k) = imresize(S0(:,:,k), 0.5);
% end
% S0 = tmp;
% image_num = 5;

% % Std 3 Training images
% image_num = 2;
% S0 = zeros(256, 256, image_num, 'single');
% S0(:,:,1) = imresize(single(rgb2gray(stdimage('monarch'))),[256,256]) / 255;
% S0(:,:,2) = imresize(single(rgb2gray(stdimage('airplane'))),[256,256]) / 255;
% % tmp = single(stdimage('man.grey')) / 255;
% % S0(:,:,3) = tmp(101:612, 101:612);
% %Reduce images size to speed up demo script


%% Initialization and Setting up Parameters
% %Filter input images and compute highpass images
npd = 16;
fltlmbd = 6;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);

% Construct initial dictionary
num_dict_high = 5; %smaller dictionary
num_dict_mid = 5;
num_dict_low = 5;
num_dict = num_dict_high+num_dict_low+num_dict_mid;
small_filter_sz = 12;
mid_filter_sz = 12;
big_filter_sz =12;

D0 = zeros(big_filter_sz,big_filter_sz,num_dict, 'single');

%initialization for high,mid,low
D0(small_filter_sz/4+1:small_filter_sz/4*3,small_filter_sz/4+1:small_filter_sz/4*3 ...
,1:num_dict_high) = single(randn(small_filter_sz/2,small_filter_sz/2,num_dict_high));

D0(mid_filter_sz/4+1:mid_filter_sz/4*3,mid_filter_sz/4+1:mid_filter_sz/4*3 ...
,num_dict_high+1:num_dict_high+num_dict_mid) = single(randn(mid_filter_sz/2,mid_filter_sz/2,num_dict_mid));

D0(big_filter_sz/4+1:big_filter_sz/4*3,big_filter_sz/4+1:big_filter_sz/4*3 ...
,num_dict_high+num_dict_mid+1:end) = single(randn(big_filter_sz/2,big_filter_sz/2,num_dict_low));

% D0(big_filter_sz/4+1:big_filter_sz/4*3,big_filter_sz/4+1:big_filter_sz/4*3 ...
% ,num_dict_high+num_dict_mid+1:end) = single(repmat(gauss2d(big_filter_sz/2,big_filter_sz/2,big_filter_sz/4),[1,1,num_dict_low]));

% Set up cbpdndliu parameters
lambda = .22;
mu1= 0; %grd on Dict
mu2 = 0;    %grd on Coef
wl1 = [ones(1,num_dict_mid),ones(1,num_dict_low)];  
wgrd = [ones(1,num_dict_mid)*.0,ones(1,num_dict_low)*.0];

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
opt.GrCoef = 0;

opt.L1Weight = reshape([ones(1,num_dict_high),wl1],1,1,num_dict);
opt.DgrdWeight = reshape([zeros(1,num_dict_high),wgrd],1,1,num_dict);
opt.XgrdWeight = opt.DgrdWeight; %remove if you don't want coefficient smoothing
opt.DictFilterSizes = [repmat([small_filter_sz;small_filter_sz],1,num_dict_high),...
    repmat([mid_filter_sz;mid_filter_sz],1,num_dict_mid),...
    repmat([big_filter_sz;big_filter_sz],1,num_dict_low)];



%% Training Dictionary
% Do dictionary learning
[D, X, optinf] = cbpdndliu_grd(D0, Sh, mu1,mu2, lambda, opt);


small_filter_sz = small_filter_sz/2;
mid_filter_sz = mid_filter_sz / 2;
big_filter_sz = big_filter_sz / 2;

D0(small_filter_sz/4+1:small_filter_sz/4*3,small_filter_sz/4+1:small_filter_sz/4*3 ...
,1:num_dict_high) = single(randn(small_filter_sz/2,small_filter_sz/2,num_dict_high));

D0(mid_filter_sz/4+1:mid_filter_sz/4*3,mid_filter_sz/4+1:mid_filter_sz/4*3 ...
,num_dict_high+1:num_dict_high+num_dict_mid) = single(randn(mid_filter_sz/2,mid_filter_sz/2,num_dict_mid));

D0(big_filter_sz/4+1:big_filter_sz/4*3,big_filter_sz/4+1:big_filter_sz/4*3 ...
,num_dict_high+num_dict_mid+1:end) = single(randn(big_filter_sz/2,big_filter_sz/2,num_dict_low));

opt.DictFilterSizes = [repmat([small_filter_sz;small_filter_sz],1,num_dict_high),...
    repmat([mid_filter_sz;mid_filter_sz],1,num_dict_mid),...
    repmat([big_filter_sz;big_filter_sz],1,num_dict_low)];

Sh = imresize(Sh,.5);
[D2,X2,optinf] = cbpdndliu_grd(D0,Sh,mu1,mu2,lambda,opt);



%% Plotting and Saving Info
% Plot functional value evolution
figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');

% Display learned dictionary
o1.grey = 1;
o1.unifscale = 1;
o1.fsz = opt.DictFilterSizes*2;
h = cell(1,4);
h{1} = dict_plot(D,o1);

o1.grey = 1;
o1.unifscale = 1;
o1.fsz = opt.DictFilterSizes;
h = cell(1,4);
h{2} = dict_plot(D2,o1);

% %plot the "high freq" plus "low freq" for all stuff
% DX = zeros(size(X,1),size(X,2),size(X,4));
% for k = 1:image_num
%     DX(:,:,k) = convsum(D,X(:,:,:,k),1:1:num_dict_high);
% end
% opt1.grey = 1;
% opt1.unifscale = 0;
% h{2} = square_plot(DX,opt1);
% 
% %mid
% DX = zeros(size(X,1),size(X,2),size(X,4));
% for k = 1:image_num
%     DX(:,:,k) = convsum(D,X(:,:,:,k),num_dict_high+1:num_dict_high+num_dict_mid);
% end
% opt1.grey = 1;
% opt1.unifscale = 0;
% h{3} = square_plot(DX,opt1);
% 
% 
% 
% %low
% DX = zeros(size(X,1),size(X,2),size(X,4));
% for k = 1:image_num
%     DX(:,:,k) = convsum(D,X(:,:,:,k),num_dict_high+num_dict_mid+1:num_dict);
% end
% opt1.grey = 1;
% opt1.unifscale = 0;
% h{4} = square_plot(DX,opt1);
% 
% 
% if ~exist(strcat(sporco_path,'/Results/DictLearn/',folder_tag),'dir')
%     mkdir(strcat(sporco_path,'/Results/DictLearn/',folder_tag));
% end
% name = {'Dict','LowRec','MidRec','HighRec'};
% for i = 1:4
%     saveas(h{i},strcat(sporco_path,'/Results/DictLearn/',folder_tag,'/',...
%         tag,name{i}),'fig');
% end
% 
% if ~strcmp(tag,'test')
%     par.mu1 = mu1;
%     par.mu2 = mu2;
%     par.lambda = lambda;
%     par.wl1 = wl1;
%     par.wgrd = wgrd;
%     par.fsize = opt.DictFilterSizes;
%     save(strcat(sporco_path,'/Results/DictLearn/',folder_tag,'/Dict',tag,'.mat'),'D','par');
% end