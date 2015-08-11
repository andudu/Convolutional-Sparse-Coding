% Testing image inpainting

% Training images 512 x 512 
%s = single(stdimage('lena.grey'))/255 ;

% Lightning image
load('lightning.mat');
s = single(x)/255;

%generate missing pixels
ind = [];
ind(1,:) = randi(size(x,1),1,ceil(size(x,1)*size(x,2)*.3));
ind(2,:) = randi(size(x,2),1,ceil(size(x,1)*size(x,2)*.3));

sn = s;
for i = 1:size(ind,2)
    sn(ind(1,i),ind(2,i)) = 1;
end




% Highpass filter test image
npd = 16;
fltlmbd = 5;
[sl, sh] = lowpass(s, fltlmbd, npd);

for i = 1:size(ind,2)
    sh(ind(1,i),ind(2,i)) = 1;
end



% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
delta = zeros(12,12);
delta(1,1) = 1;
D(:,:,end+1) = delta;
numdict = size(D,3);

opt.L1Weight = ones(size(x,1),size(x,2),37);
opt.L1Weight(:,:,37) = 100*ones(size(x));
for i = 1:size(ind,2)
    opt.L1Weight(ind(1,i),ind(2,i),37) = 0;
end

lambda = .03;
opt.Verbose = 1;
opt.MaxMainIter = 250;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.AutoRhoPeriod = 4;


[X,~] = cbpdn(D,sh,lambda,opt);
sh_rec = convsum(D,X,1:1:numdict-1);


