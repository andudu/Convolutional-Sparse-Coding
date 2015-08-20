% small scale denoise test 

%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting options for Graph Laplacian
optl = {};
optl.wsz = [60,60];
optl.psz = [8,8];
optl.neig = 40;
optl.Lformat = 'Sparse';
optl.Laplacian = 'u';
optl.Graph.tau = 1;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'Window';
optl.Graph.nsz = [8,8];
optl.Graph.k = [];


imsz = optl.wsz+optl.psz-[1,1] ;
psz = optl.psz;
stpsz = [1,1];
load('stdnoise.mat');
n = r_noise(1:imsz(1),1:imsz(2));
Dversion = 'lenapatch+noise';
Wversion = 0;
%%%%%%%%%%%%%%%%%%   Generating Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(Dversion, 'lenapatch'),
    %patch from lena
    s = double(stdimage('lena.grey')) / 255;
    s = imresize(s,.5);
    s_ref = s(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
    sn = s(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
    [sl,sh] = lowpass(s,7,15);
    s = sh(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
    sl = sl(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
end

if strcmp(Dversion, 'lenapatch+noise'),
    %patch from lena
    s = double(stdimage('lena.grey')) / 255;
    s = imresize(s,.5);
    s_ref = s(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
    s = s+r_noise;
    sn = s(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
    [sl,sh] = lowpass(s,7,15);
    s = sh(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
    sl = sl(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
end
scrop = s(1:wsz(1),1:wsz(2));
    


%%%%%%%%%%%%%%%%%   Building the Graph Weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1d = imageblocks(s,optl.psz,[1,1]);
s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));
optl.Graph.coefsz = optl.wsz;
W = graphgen(s1d, optl.Graph);
if optl.Laplacian == 'n'
    L = nlap(W);
else
    L = ulap(W);
end

%%%%%%%%%%%%%%%%%%%%%%%%%   Generating L1 Weights %%%%%%%%%%%%%%%%%%%%%%%%%



% % use some sort of diffusion trick. Not well thought up at the moment
% ind = [26,49]-[6,6];
% flatind = rec2flat(ind,imsz,psz,stpsz);
% 
% Ltemp = L;
% dt = .5;
% x = zeros(size(Ltemp,1),1);
% x(flatind) = .1;%initial condition
% 
% for k = 1:6
% y = x - dt*Ltemp*x;
% wsz = optl.wsz;
% x = y;
% y = reshape(y,wsz(1),wsz(2))';
% figure;
% imagesc(y);
% colorbar;
% end;


% use a degree weighting scheme on the cbpdn coefficient

mu = .1;
lambda = 0.3;
opt = {};
opt.Verbose = 0;
opt.MaxMainIter = 50;
opt.rho = 10;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
load('CacheData/Dict_12x12.mat');
[Xcn,~] = cbpdngr(D,scrop,lambda,mu,opt);

d = sum(abs(W),2);
d = d./max(d);
wsz = optl.wsz;
d = reshape(d,wsz(1),wsz(2))';

X = bsxfun(@times,abs(Xcn),d);
for i = 1:size(Weight,3)
    temp = abs(X(:,:,i));
    d = max(vec(temp));
    temp = temp/d;
    Weight(:,:,i) = exp(-1.8*temp);      
end
square_plot(Weight,{});


%%%%%%%%%%%%%%%%%%%%%%%%%   Testing the Weighting %%%%%%%%%%%%%%%%%%%%%%%%%

lambda = .4;
opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 100;
opt.rho = 10;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = Weight;

[X1,~] = cbpdn(D,scrop,lambda,opt);

lambda = .3;
opt.L1Weight = 1;
[X2,~] = cbpdn(D,scrop,lambda,opt);

% reconstructing conv
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

                           
                           
sl = sl(1:wsz(1), 1:wsz(2));    
s_ref = s_ref(1:wsz(1), 1:wsz(2));  
shrec1 = scnv(D,X1);
srec1 = shrec1 + sl;
p1 = psnr(srec1,s_ref);
disp(['psnr from Reweigh: ',num2str(p1)]);

figure; 
imagesc(srec1);
title(['Reweighed psnr = ',num2str(p1)]);
colormap(gray);
colorbar;

figure; 
imagesc(shrec1);
title(['Reweighed sh']);
colormap(gray);
colorbar;


figure; 
imagesc(sum3(X1));
title(['Reweighed coeff']);
colormap(gray);
colorbar;



shrec1 = scnv(D,X2);
srec1 = shrec1 + sl;
p1 = psnr(srec1,s_ref);
disp(['psnr from cbpdn: ',num2str(p1)]);

figure; 
imagesc(srec1);
title(['cbpdn psnr = ',num2str(p1)]);
colormap(gray);
colorbar;

figure; 
imagesc(shrec1);
title(['cbpdn sh']);
colormap(gray);
colorbar;


figure; 
imagesc(sum3(X2));
title(['cbpdn coeff']);
colormap(gray);
colorbar;




















