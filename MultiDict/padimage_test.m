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
[~,S0] = lowpass(S0, fltlmbd, npd);

%Reduce images size to speed up demo script
S1 = zeros(512,512,5);
for k = 1:size(S0,3),
  a = imresize(S0(:,:,k), [256+7,256+7]); % pad a little bit to reduce edge effects
  [~,a] = lowpass(a, fltlmbd, npd);
  S1(:,:, k) = padarray(a,[512-size(a,1),512-size(a,2)],'post'); 
end
% Filter input images and compute highpass images
npd = 16;
fltlmbd = 5;
[~, Sh] = lowpass(S0, fltlmbd, npd);


%Reduce images size to speed up demo script
S2 = zeros(512, 512, 5, 'single');
for k = 1:size(S0,3),
  a = imresize(S0(:,:,k), [128+7,128+7]); % pad a little bit to reduce edge effects
  [~,a] = lowpass(a, fltlmbd, npd);
  S2(:,:, k) = padarray(a,[512-size(a,1),512-size(a,2)],'post'); 
end

S = zeros(size(S0,1),size(S0,2),3*size(S0,3));
opt.imsz = [];
for k = 1:size(S0,3)
    S(:,:,k) = S0(:,:,k);
    opt.imsz = [opt.imsz, [512;512]];
end
for k = 1:size(S0,3)
    S(:,:,k+size(S0,3)) = S1(:,:,k);
    opt.imsz = [opt.imsz, [256;256]];
end    
for k = 1:size(S0,3)
    S(:,:,k+2*size(S0,3)) = S2(:,:,k);
    opt.imsz = [opt.imsz, [128;128]];
end 



% % Construct initial dictionary
D0 = zeros(8,8,15, 'single');
D0(3:6,3:6,:) = single(randn(4,4,15));

% Load a standard dictionary
% load([sporco_path '/Data/ConvDict.mat']);
% dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
% f = dmap('8x8x32');
% D0 = f(:,:,1:12);

% Set up cbpdndliu parameters
lambda = 0.2;
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(S,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.2;
opt.DRelaxParam = 1.2;
opt.StdResiduals = 0;

% Do dictionary learning
[D, X, optinf] = padimage_multilearn(D0, S, lambda, opt);

tag = '5im15dictcold3Layer';
% Display learned dictionary
o1.grey =1;
o1.unifscale =0;
square_plot(D,o1);
saveas(gcf,['PadResults/',tag,'_dict'],'fig');
o1.grey =0;
a = reshape(sum(abs(X),3),size(X,1),size(X,2),size(X,4));
square_plot(a,o1);
saveas(gcf,['PadResults/',tag,'_total'],'fig');
save(['PadResults/',tag,'_dict.mat'],'D');


% Plot functional value evolution
figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');