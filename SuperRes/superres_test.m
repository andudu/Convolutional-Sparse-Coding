%test script for super resolution


% Construct lowpass filtering and downsampling operators
g = gauss2d([5 5], 0.75);
smooth = @(x) imfilter(x, g, 'symmetric');
dnsmpl = @(x) x(1:2:end, 1:2:end);
dwnres = @(x) dnsmpl(smooth(x));

%Training on flicker images
flicker_ind = [4,5,7,10,11,12];
load('Flicker1_512_split.mat');
S0 = [];
k = 0;
for i = flicker_ind
    S0(:,:,(k)*4+1) = S(:,:,(i-1)*4+1);
    S0(:,:,(k)*4+2) = S(:,:,(i-1)*4+2);    
    S0(:,:,(k)*4+3) = S(:,:,(i-1)*4+3);
    S0(:,:,(k)*4+4) = S(:,:,(i-1)*4+4);
    k = k+1;
end
S0 = single(S0)/255;
image_num = length(flicker_ind);
clear S;

%Downsample image
S1 = zeros(128, 128, 4*image_num, 'single');
for k = 1:size(S0,3),
    S1(:,:,k) = dwnres(S0(:,:,k));
end

% Filter input images and compute highpass images
npd = 16;
fltlmbd = 5;
[Sl, Sh] = lowpass(S0, fltlmbd, npd);
[Sl1, Sh1] = lowpass(S1, fltlmbd, npd);

% % Construct initial dictionary
numdict = 18;
D0 = zeros(12,12,numdict, 'single');
D0(4:9,4:9,:) = single(randn(6,6,numdict));
D0_bar = zeros(6,6,numdict, 'single');
D0_bar(2:5,2:5,:) = single(randn(4,4,numdict));

% load([sporco_path '/Data/ConvDict.mat']);
% dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
% f = dmap('8x8x32');
% D0 = f(:,:,1:12);

% Set up cbpdndliu parameters
lambda = 0.2;
opt = [];
opt.Verbose = 1;
opt.rho = 50*lambda + 0.5;
opt.sigma = size(Sh,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.AutoDelta = 1;
opt.beta = 1.4;

%just add a seperate stepsize for all _bar variable
opt.AutoRho_bar = 1;
opt.AutoRhoPeriod_bar = 10;
opt.AutoSigma_bar = 1;
opt.AutoSigmaPeriod_bar = 10;
opt.AutoDelta = 1;
opt.AutoDelta = 5;
opt.delta = 3;
opt.XRelaxParam = 1.5;
opt.DRelaxParam = 1.5;
opt.LinSolve = 'SM';

% Do dictionary learning
opt.MaxMainIter = 30;
[D, D_bar,X,X_bar, optinf] = superres(D0, D0_bar, Sh, Sh1, lambda, opt);

tag = ['Flicker',num2str(image_num),'im',num2str(numdict),'dict','cold'];
% Display learned dictionary
o1.grey =1;
o1.unifscale =0;
square_plot(D,o1);
saveas(gcf,['SuperresResults/',tag,'_highresdict'],'fig');

% Display learned dictionary
o1.grey =1;
o1.unifscale =0;
square_plot(D_bar,o1);
saveas(gcf,['SuperresResults/',tag,'_lowresdict'],'fig');
o1.grey =0;
a = reshape(sum(abs(X),3),size(X,1),size(X,2),size(X,4));
square_plot(a,o1);

saveas(gcf,['SuperresResults',tag,'_coeff'],'fig');

square_plot(Sh,o1);

square_plot(Sh1,o1);

save(['SuperresResults',tag,'_dict.mat'],'D','D_bar','X');

% num = 2;
% a = reshape(sum(abs(X(:,:,num,:)),3),size(X,1),size(X,2),size(X,4));
% square_plot(a,o1);
% a = reshape(sum(abs(X_bar(:,:,num,:)),3),size(X_bar,1),size(X_bar,2),size(X_bar,4));
% square_plot(a,o1);


% for i = 1:5
%     temp = convsum(D,X(:,:,:,i),1:14);
%     figure;
%     imagesc(temp);
%     temp = convsum(D_bar,X(:,:,:,i),1:14);
%     figure;
%     imagesc(temp);
%     figure;
%     imagesc(Sh(:,:,i));
%     
% end


% Plot functional value evolution
figure;
plot(optinf.itstat(:,2));
xlabel('Iterations');
ylabel('Functional value');