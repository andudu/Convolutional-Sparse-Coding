%code for testing for patch based denoising

%parameters
patchSize = 12;
stepSize = 1;
imsz = 256;
tag = '2noise+low1e-1lam2e-1_';


%%%%%%%%%%%%%%%%%%%%%%%%% data generation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pure noise + big line + little line
load ('noise_data.mat');
s_ref = zeros(size(s));
% s_ref(120:121,:) = 1;  %big line
% s_ref(150:151,120:124) = 1;  %small line
% s_ref(90:91,:) = .5;  %big line with half intensity
% s_ref(150:151,140:144) = .5; %small line with half intensity
% s_ref(150:151,100:104) = .2; %small line with half intensity
s = s+s_ref;
[sl,sh] = lowpass(s,10,16);



%noise + lena;

% s_ref = single(rgbtogrey(stdimage('lena')))/255;
% s_ref = imresize(s_ref,.5);
% s = s+s_ref;
% [sl,sh] = lowpass(s,5,16);

%%%%%%%%%%%%%%%%%%%%%%% preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s_flat = imageblocks(s, [patchSize patchSize]);
% s_flat = reshape(s_flat, size(s_flat,1)*size(s_flat,2), size(s_flat,3));
% s_flat_mean = mean(s_flat,1);
% s1d = bsxfun(@minus, s_flat, s_flat_mean);


% temp = repmat(s_flat_mean,[patchSize*patchSize,1]);
% s_flat_2d = unpatchify(s_rec, patchSize, stepSize, [imsz,imsz]);


%take it directly from the lowpass
s_flat = imageblocks(sh, [patchSize patchSize]);
s_flat = reshape(s_flat, size(s_flat,1)*size(s_flat,2), size(s_flat,3));
s1d = s_flat;


% Construct initial dictionary
% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D1d = reshape(D,size(D,1)*size(D,2),size(D,3));


% Set up bpdndliu parameters
lambda = .2;
%lambda_patch = lambda;
opt = [];
opt.Verbose = 1;
opt.rho = 100;
opt.MaxMainIter =200;
opt.RelStopTol = 1e-6;




%%%%%%%%%%%%%%%%%%%%%%%%%%%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running optimization
% [x1, ~] = bpdn(D1d, s1d, lambda, opt);

opt.MaxMainIter =150;
opt.RelStopTol = 1e-3;
[x2, ~] = cbpdn(D, sh, lambda, opt);

scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');
shconv = scnv(D,x2);
scconv = shconv+sl;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstructing patch
% s_rec1 = D1d*x1;
% s_rec = bsxfun(@plus,s_rec1,s_flat_mean);
% srec_patch = unpatchify(s_rec, patchSize, stepSize, [imsz,imsz]);
% srec_high = unpatchify(s_rec1, patchSize, stepSize, [imsz,imsz]);


scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');
shconv = scnv(D,x2);
scconv = shconv+sl;



%plot the reconstruction from each element
temp = zeros(size(x2));
for i = 1:36
    temp(:,:,i) = convsum(D,x2,i);
end

o.unifscale = 1;
square_plot(temp,o);


% figure;
% imagesc(srec_patch); colorbar;
% title('srec_patch');
%saveas(gcf,'patchrec','fig');

figure;
imagesc(scconv);
title('recconv');
colorbar;
saveas(gcf,'convrec','fig');


% 
% %reconstructing coefficient maps
% x1_2d = zeros(imsz,imsz,size(D,3));
% for k = 1:size(x1,1)
%     k
%     for i = 1:size(x1,2)
%         ind = flat2rec(i,imsz,patchSize);
%         x1_2d(ind(1),ind(2),k) = x1(k,i);
%     end
% end
%     
% 
% 
% 
% %reconstruction of the image from coeff 6
% foo = D1d(:,6) * x1(6,:);
% s4 = unpatchify(foo,patchSize,stepSize,[imsz,imsz]);
% 
% 
% 
% 




