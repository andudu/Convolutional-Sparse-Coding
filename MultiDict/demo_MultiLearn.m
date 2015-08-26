%Testing the Multilearn Dictionary Stuff. 


% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D0 = dmap('8x8x32');
D0 = D0(:,:,1:15);
num_dict = size(D,3);

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
fltlmbd = 5;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);




% Set up cbpdndliu parameters
lambda = 0.23;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 250;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(Sh,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 6;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 6;
opt.XRelaxParam = 1.6;
opt.DRelaxParam = 1.6;
opt.ImSize = [512;512]*[.5];


% Do dictionary learning
[D, X, optinf] = cbpdndliu(D0, Sh1, lambda, opt);


% Display learned dictionary
figure;
imdisp(tiledict(D));

% Plot functional value evolution
figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');

sx = sum(abs(X),3);
x1 = sx(:,:,1,1);
figure;
imagesc(x1);

