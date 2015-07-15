%this script tests random zero fillings 

% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('8x8x32');
num_dict = size(D,3);

%loda a test image
imagename = 'lena';
s = single(rgbtogrey(stdimage(imagename)))/255; %use barbara or lena here
s = imresize(s, 0.5);
s_ref = s;
%add the noise (10 percent peak signal)
sigma = 20/255;
s = s+ sigma*randn(size(s));



%function for convsum
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

% Highpass filter test image
npd = 16;
fltlmbd = 5;
[sl, sh] = lowpass(s, fltlmbd, npd);

%setting up params
lambda = 0.13; %best lambda
opt = [];
opt.Verbose = 0;
opt.MaxMainIter = 200;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;


% %compute sparse representation
% disp('Regular CBPDN');
% [X_d, ~] = cbpdn_zfill(D, sh, lambda,z_ind, opt);
% sh_d = scnv(D,X_d);
% s_d = sh_d + sl;
% ps_nr_d = psnr(s_d,s_ref);
% figure;
% imagesc(s_d); colormap(gray);
% title(['psnr of Regular = ',num2str(ps_nr_d)]);


%do an averaging
ratio = .1;
a = zeros(size(s));
for i = 1:20
    disp(num2str(i));
    z_ind = randi(256,[2,ceil(ratio*256*256)]);
    [X_zfill, ~] = cbpdn_zfill(D, sh, lambda,z_ind, opt);
    sh_zfill = scnv(D,X_zfill);
    s_zfill = sh_zfill + sl;
    a = a+s_zfill;
end
a = a/20;
figure;
imagesc(a); colormap(gray);
title(['psnr of averaged = ',num2str(psnr(a,s_ref))]);







