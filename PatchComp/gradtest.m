%code for testing Averaging Effect of a spatial gradient term. 


%%%%%%%%%%%%%%%%%%%%%%%%% data generation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %pure noise 
% s = .1*randn(imsz,imsz); 
% s_ref = zeros(size(s));


% %pure noise 
% sn = .1*randn(imsz,imsz); 
% s_ref = zeros(size(sn));
% s = s_ref+sn;


%lena
load('noise_data.mat');
s_ref = single(rgbtogrey(stdimage('lena')))/255;
s_ref = imresize(s_ref,.5);
s = s_ref+s;
imsz = size(s,1);


%%%%%%%%%%%%%%%%%%%%%%% preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sl,sh] = lowpass(s,5,16);

% Construct initial dictionary
% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');


% Set up bpdndliu parameters
lambda = .2;
opt = [];
opt.Verbose = 1;
opt.rho = 100;
opt.MaxMainIter =150;
opt.RelStopTol = 1e-3;
mu = 5;



%%%%%%%%%%%%%%%%%%%%%%%%%%%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running optimization
[x1, ~] = cbpdn(D, sh, lambda, opt);
lambda = .01;
[x2,~] = cbpdngr(D,sh,lambda,mu, opt);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstructing conv
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');
shrec1 = scnv(D,x1);
shrec2 = scnv(D,x2);

figure;
imagesc(sum3(abs(x1)));
colorbar;
title('Coeff l1');

figure;
imagesc(sum3(abs(x2)));
colorbar;
title('Coeff grad');

figure;
imagesc(sum3(abs(shrec1)));
colorbar;
title('Highrec l1');

figure;
imagesc(sum3(abs(shrec2)));
colorbar;
title('Highrec grad');


figure;
imagesc(sum3(abs(shrec1+sl)));
colorbar;
title('Rec l1');

figure;
imagesc(sum3(abs(shrec2+sl)));
colorbar;
title('Rec grad');
















