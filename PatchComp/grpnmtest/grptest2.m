% script for testing group norm

% Load data
% s = zeros(256,256);
% box_lu = [124,124];
% box_rd = [136,136];
% line_l = 130;
% line_r = 137;
% line_h = 135;
% line_width = 2;
% 
% s(line_h:line_h+line_width-1,line_l:line_r) = 1; %short line
% 
% s(100:101,:) = 1;

% load('noisepatch.mat');
% s = zeros(256,256);
% s(100:120,100:115) = a;

load('noise_data.mat');


% Construct initial dictionary
% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D = D(:,:,6);

%[sl,sh] = lowpass(s,5,15);
sh = s;


opt = [];
opt.Verbose = 1;
opt.rho = 150;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 5;
opt.MaxMainIter =100;
opt.RelStopTol = 1e-3;
opt.L1Weight = ones(size(sh));
opt.L1Weight(236:242,141:145) = 0;

lambda = .2;
mu = .2;
g = zeros(size(sh));
g(236:242,141:145)= 1; %group the box into one 


[x1,~] = cbpdngrp(D,sh,lambda,mu,g,opt);
opt.L1Weight = 1;
[x2,~] = cbpdn(D,sh,lambda,opt);


% Visualization
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');
shrec1 = scnv(D,x1);
shrec2 = scnv(D,x2);
srec1 = shrec1;
srec2 = shrec2;

figure;
imagesc(srec1);
title('reconstruction from group norm');
colorbar;
saveas(gcf,'gprec','fig');

figure; 
imagesc(sum3(abs(x1)));
title('sum of coeff gp norm');
colorbar;
saveas(gcf,'gpcoef','fig');

figure;
imagesc(srec2);
title('reconstruction from L1 norm');
colorbar;

saveas(gcf,'l1rec','fig');

figure; 
imagesc(sum3(abs(x2)));
title('sum of coeff L1 norm');
colorbar;

saveas(gcf,'l1coef','fig');






