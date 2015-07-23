% script for testing group norm

% Load data
s_ref = single(rgbtogrey(stdimage('lena')))/255;
s_low = imresize(s_ref,.5);


% Construct initial dictionary
% Load dictionary
%load('Flicker30im34dictcold1_dict.mat');
load('Flicker20im34dictcold_dict.mat');

[sl_low,sh_low] = lowpass(s_low,5,15);
[sl_ref, sh_ref] = lowpass(s_ref,5,15);

opt = [];
opt.Verbose = 1;
opt.rho = 150;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 3;
opt.MaxMainIter =100;
opt.RelStopTol = 1e-3;
lambda = .01;


[x2,~] = cbpdn(D_bar,sh_low,lambda,opt);


% Visualization
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

x = PxT(x2); 


                           
shrec = scnv(D,x);

figure;
imagesc(shrec);
title('superres high');
colorbar;

figure;
imagesc(sh_ref);
title('Original high');
colorbar;


[~, foo] = lowpass(imresize(s_low,2),4,15);
figure;
imagesc(foo);
title('imresize');
colorbar;

foo = zeros(size(x));
for i = 1:34
   foo(:,:,i) = convsum(D,x,i); 
end
square_plot(foo,{});
title('component rec');


    


