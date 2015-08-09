% Testing image inpainting

% Training images 512 x 512 
s = single(stdimage('lena.grey'))/255 ;

%generate missing pixels
ind = randi(512,2,ceil(512*512*.3));

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

opt.L1Weight = ones(512,512,37);
opt.L1Weight(:,:,37) = 100*ones(512,512);
for i = 1:size(ind,2)
    opt.L1Weight(ind(1,i),ind(2,i),37) = 0;
end

lambda = .03;
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.AutoRhoPeriod = 4;


[X,~] = cbpdn(D,sh,lambda,opt);
sh_rec = convsum(D,X,1:1:numdict-1);

snoise = convsum(D,X,numdict);


% 
% figure;
% imagesc(sn);
% axis off;
% saveas(gcf,'saltpeppernoisylena','png');
% 
% figure;
% imagesc(sh_rec+sl);
% axis off;
% saveas(gcf,'saltpepperreclena','png');

imwrite(sn,'iplenaoriginal.png');
imwrite(sh_rec+sl,'iplenarec.png');
