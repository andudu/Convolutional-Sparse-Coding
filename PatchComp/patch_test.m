%code for testing for patch based denoising

%parameters
patchSize = 12;
stepSize = 1;
imsz = 256;
tag = '2noise+low1e-1lam2e-1_';


%%%%%%%%%%%%%%%%%%%%%%%%% data generation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %pure noise 
% s = .1*randn(imsz,imsz); 
% s_ref = zeros(size(s));


% %pure noise 
% sn = .1*randn(imsz,imsz); 
% s_ref = zeros(size(sn));
% s = s_ref+sn;

%save('noise_data','s');

%lena
s = single(rgbtogrey(stdimage('lena')))/255;
s = imresize(s,.5);
s_ref = s;
s = s+.1*randn(size(s));
imsz = size(s,1);


%%%%%%%%%%%%%%%%%%%%%%% preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_flat = imageblocks(s, [patchSize patchSize]);
s_flat = reshape(s_flat, size(s_flat,1)*size(s_flat,2), size(s_flat,3));
s_flat_mean = mean(s_flat,1);
s1d = bsxfun(@minus, s_flat, s_flat_mean);
[sl,sh] = lowpass(s,5,16);

a = zeros(imsz,imsz);
for i = 1:size(s_flat,2)
   ind = flat2rec(i,imsz,patchSize); 
   a(ind(1):ind(1)+patchSize-1,ind(2):ind(2)+patchSize-1) =a(ind(1):ind(1)+patchSize-1,ind(2):ind(2)+patchSize-1) + reshape(s_flat(:,i),patchSize,patchSize);
end



% Construct initial dictionary
% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D1d = reshape(D,size(D,1)*size(D,2),size(D,3));


% Set up bpdndliu parameters
lambda = .5;
lambda_patch = (.3/2)*lambda;
%lambda_patch = lambda;
opt = [];
opt.Verbose = 1;
opt.rho = 100;
opt.MaxMainIter =200;
opt.RelStopTol = 1e-6;




%%%%%%%%%%%%%%%%%%%%%%%%%%%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running optimization
[x1, ~] = bpdn(D1d, s1d, lambda_patch, opt);
opt.MaxMainIter =150;
opt.RelStopTol = 1e-3;
[x2, ~] = cbpdn(D, sh, lambda, opt);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstructing patch
s_rec1 = D1d*x1;
s_rec = bsxfun(@plus,s_rec1,s_flat_mean);
srec_patch = unpatchify(s_rec, patchSize, stepSize, [imsz,imsz]);
srec_high = unpatchify(s_rec1, patchSize, stepSize, [imsz,imsz]);

% reconstructing conv
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');
srec_conv = scnv(D,x2)+sl;

% plotting
figure;
subplot(1,2,1);
imagesc(s_ref);
title('Original');
colorbar;
subplot(1,2,2);
imagesc(s);
title('Noisy');
colorbar;
saveas(gcf,['DenoiseResult/',tag,'rawim'],'fig');


figure;
subplot(1,2,1);
imagesc(srec_high);
title('Highfreq-patch');
colorbar;
subplot(1,2,2);
imagesc(srec_conv-sl);
title('Highfreq-conv');
colorbar;
saveas(gcf,['DenoiseResult/',tag,'Highfreq'],'fig');



figure;
subplot(1,2,1);
imagesc(srec_patch-srec_high);
title('Lowfreq-patch');
colorbar;
subplot(1,2,2);
imagesc(sl);
title('Lowfreq-conv');
colorbar;
saveas(gcf,['DenoiseResult/',tag,'Lowfreq'],'fig');



figure;
subplot(1,2,1);
imagesc(srec_patch);
title('Recon-patch');
colorbar;
subplot(1,2,2);
imagesc(srec_conv);
title('Recon-conv');
colorbar;
saveas(gcf,['DenoiseResult/',tag,'Recon'],'fig');

figure;
subplot(1,2,1);
imagesc(srec_patch-s_ref);
title('Diff-patch');
colorbar;
subplot(1,2,2);
imagesc(srec_conv-s_ref);
title('Diff-conv');
colorbar;
saveas(gcf,['DenoiseResult/',tag,'Diff'],'fig');


