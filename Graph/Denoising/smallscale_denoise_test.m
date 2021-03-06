% small scale denoise test 

%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting options for Graph Laplacian
optl = {};
optl.wsz = [60,60];
optl.psz = [8,8];
optl.neig = 40;
optl.Lformat = 'Eig';
optl.Laplacian = 'n';
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
Dversion = 'simpleline+noise';
Wversion = 0;
%%%%%%%%%%%%%%%%%%   Generating Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Dversion, 'simpleline'),
    % simple line
    s = zeros(imsz);
    s(35:36,:) = 1;
    s_ref = s;
    [sl,sh] = lowpass(s,7,15);
    sn = s;
    s = sh;
end

if strcmp(Dversion, 'simpleline+noise'),
    % simple line
    s = zeros(imsz);
    s(35:36,:) = 1;
    s_ref = s;
    s = s+n;
    [sl,sh] = lowpass(s,7,15);   
    sn = s;
    s = sh;
end

if strcmp(Dversion, 'longshortline'),
     % simple line
    s = zeros(imsz);
    s(35:36,:) = 1;
    s(45:45,32:36) = 1;
    s_ref = s;
    [sl,sh] = lowpass(s,7,15);
    sn = s;
    s = sh;
end


if strcmp(Dversion, 'longshortline+noise'),
     % simple line
    s = zeros(imsz);
    s(35:36,:) = 1;
    s(45:45,32:36) = 1;
    s_ref = s;
    s = s+n;
    [sl,sh] = lowpass(s,7,15);    
    sn = s;
    s = sh;
end


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

    


%%%%%%%%%%%%%%%%%   Building the Graph Weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version 1. Connect nonzero/zero elements to its nonzero/zero neighbors
% Perfect paring in this case
if Wversion == 1
    W = sparse(optl.wsz(1)*optl.wsz(2), optl.wsz(1)*optl.wsz(2));
    for i = 1:optl.wsz(1)
        for j = 1:optl.wsz(2)
            flind = rec2flat([i,j],imsz,psz,stpsz);
            if (any(s_ref(i:i+psz(1)-1,j:j+psz(2)-1) == 1)), % actually a patch containing a line
                for k2 = -5:5 %connect to its 5 horizontal neighbors
                    if  j+k2>0  && j+k2<= optl.wsz(2)
                        flind2 = rec2flat([i,j+k2],imsz,psz,stpsz);
                        W(flind,flind2) = 1;
                    end
                end
            else
                for k1 = -1:1
                    for k2 = -1:1
                        if i+k1>0 && j+k2>0 && i+k1<=optl.wsz(1) && j+k2<= optl.wsz(2)
                            if ~any(s_ref(i+k1:i+k1+psz(1)-1,j+k2:j+k2+psz(2)-1) == 1), %connect only to 0 patches
                                flind2 = rec2flat([i+k1,j+k2],imsz,psz,stpsz);
                                W(flind,flind2) = 1;
                            end
                        end
                    end
                end
            end
        end
    end
    W = max(W,W');
    if strcmp(optl.Laplacian,'n')
        L = nlap(W);
    else
        L = ulap(W);
    end
    ind1 = [1,1];
    ind2 = optl.wsz;
    temp = L;
    L = {}; L{1}.ind1 = ind1;
    L{1}.ind2 = ind2; L{1}.M = temp;
    L_perfect = L;    
end


% Version 0 Choose a scheme specified by optl
% 
% Wversion = 0;
% if Wversion == 0
%     [L,scrop] = laplacian_from_image(s,optl);
% end
% 
% 
% E1 = []; E2 = [];
% for i = 1:optl.neig
%     x = (i/optl.neig);
%     E1(i) = x^3;
%     E2(i) = 1-(x-1)^3;
% end
% 
% L1 = L;
% L2 = L;
% 
% L1{1}.E = E1';
% L2{1}.E = E2';

% %%%%%%%%%%%%%%%%%%  Testing by Denoising Experiment %%%%%%%%%%%%%%%%%%%%%%%
% 
% % cbpdn with graph regularization
mu = .05;
lambda = .15;

opt = {};
opt.Verbose = 1;
opt.MaxMainIter = 200;
opt.rho = 100*lambda + 1;
opt.sigma = opt.rho;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;
opt.L1Weight = 1;



% Load dictionary
load('CacheData/Dict_12x12.mat');
D = double(D);

opt.Lformat = optl.Lformat;
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');





% 
% % reconstructing conv
% figure; imagesc(sn); colorbar; title('noisy image'); colormap(gray);
% figure; imagesc(s); colorbar; title('noisy image high'); colormap(gray);



% %nonlocal
% [Xnl,~]= cbpdnl_lasso(D,s,L,lambda,mu,opt);
% shrecnl = scnv(D,Xnl);
% % figure; imagesc(shrecnl); colorbar; title('nonlocal high rec'); colormap(gray);
% % figure; imagesc(shrecnl+sl); colorbar; title('nonlocal rec'); colormap(gray);
% figure; imagesc(sum3(abs(Xnl))); colorbar; title('nonlocal coef');colormap(gray);
% % % 
% 
% 
% %nonlocal
% [Xnl,~]= cbpdnl_lasso(D,s,L2,lambda,mu,opt);
% shrecnl = scnv(D,Xnl);
% % figure; imagesc(shrecnl); colorbar; title('nonlocal high rec'); colormap(gray);
% % figure; imagesc(shrecnl+sl); colorbar; title('nonlocal rec'); colormap(gray);
% figure; imagesc(sum3(abs(Xnl))); colorbar; title('nonlocal coef');colormap(gray);
% % % 


% % % regular
% [Xcn,~] = cbpdn(D,s,lambda,opt);
% shrec = scnv(D,Xcn);
% figure; imagesc(shrec); colorbar; title('conv high rec'); colormap(gray);
% figure; imagesc(shrec+sl); colorbar; title('conv rec'); colormap(gray);
% figure; imagesc(sum3(abs(Xcn))); colorbar; title('conv coef');colormap(gray);
% % 


% % 
% % 
% % % nonlocal perfect
% % [Xnlp,~]= cbpdnl_lasso(D,s,L_perfect,lambda,mu,opt);
% % shrecnl = scnv(D,Xnlp);
% % figure; imagesc(shrecnl); colorbar; title('nonlocal perfect high rec'); colormap(gray);
% % figure; imagesc(shrecnl+sl); colorbar; title('nonlocal perfect rec'); colormap(gray);
% % figure; imagesc(sum3(abs(Xnlp))); colorbar; title('nonlocal perfect coef');colormap(gray);
% % % % 
% % 
% % 
% % 
% 
% % % gradient
% [Xgrd,~] = cbpdngr(D,s,lambda,mu,opt);
% shrec = scnv(D,Xgrd);
% % figure; imagesc(shrec); colorbar; title('conv grd high rec'); colormap(gray);
% % figure; imagesc(shrec+sl); colorbar; title('conv grd rec'); colormap(gray);
% figure; imagesc(sum3(abs(Xgrd))); colorbar; title('conv grd coef');colormap(gray);
% % % 
% % 
% 
% % TV regularizer
% opt.AutoSigma = 1;
% opt.AutoRho = 1;
% [Xtv,~] = cbpdntv(D,s,lambda,mu,opt);
% shrec = scnv(D,Xtv);
% % figure; imagesc(shrec); colorbar; title('conv TV high rec'); colormap(gray);
% % figure; imagesc(shrec+sl); colorbar; title('conv TV rec'); colormap(gray);
% figure; imagesc(sum3(abs(Xtv))); colorbar; title('conv TV coef');colormap(gray);


% % x2
% opt.HighMemSolve = 1;
% [Xx2,~] = cbpdnx2(D,s,lambda,mu,opt);
% shrec = scnv(D,Xx2);
% % figure; imagesc(shrec); colorbar; title('conv x2 high rec'); colormap(gray);
% % figure; imagesc(shrec+sl); colorbar; title('conv x2 rec'); colormap(gray);
% figure; imagesc(sum3(abs(Xx2))); colorbar; title('conv x2 coef');colormap(gray);
% % % % 
% 

% Reweighted L1
opt.L1Weight = 1.5*ones(size(s,1),size(s,2),size(D,3));
opt.L1Weight(30:32,19,:) = 1;
lambda = lambda + .05;
[Xcn,~] = cbpdn(D,s,lambda,opt);
lambda = lambda - .05;
shrec = scnv(D,Xcn);
figure; imagesc(shrec); colorbar; title('conv reweighl1 high rec'); colormap(gray);
figure; imagesc(shrec+sl); colorbar; title('conv reweighl1 rec'); colormap(gray);
figure; imagesc(sum3(abs(Xcn))); colorbar; title('conv reweighl1 coef');colormap(gray);
% 


% % % gradient + Reweighted L1
% opt.L1Weight = 1.5*ones(size(s,1),size(s,2),size(D,3));
% opt.L1Weight(30:32,19,:) = 1;
% [Xgrd,~] = cbpdngr(D,s,lambda,mu,opt);
% shrec = scnv(D,Xgrd);
% figure; imagesc(shrec); colorbar; title('conv grd l1 high rec'); colormap(gray);
% figure; imagesc(shrec+sl); colorbar; title('conv grd rec'); colormap(gray);
% figure; imagesc(sum3(abs(Xgrd))); colorbar; title('conv grd l1 coef');colormap(gray);

% % TV regularizer + Reweighed L1
opt.AutoSigma = 1;
opt.AutoRho = 1;
opt.L1Weight = 1.5*ones(size(s,1),size(s,2),size(D,3));
opt.L1Weight(30:32,19,:) = 1;
[Xtv,~] = cbpdntv(D,s,lambda,mu,opt);
shrec = scnv(D,Xtv);
figure; imagesc(shrec); colorbar; title('conv TV high rec'); colormap(gray);
figure; imagesc(shrec+sl); colorbar; title('conv TV rec'); colormap(gray);
figure; imagesc(sum3(abs(Xtv))); colorbar; title('conv TV Reweigh coef');colormap(gray);










