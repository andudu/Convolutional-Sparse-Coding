%script for testing simple graphs on an image of a line.

%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting options for Graph Laplacian
optl = {};
optl.wsz = [60,60];
optl.psz = [12,12];
optl.neig = 50;
optl.Lformat = 'Full';
optl.Laplacian = 'n';
optl.Graph.tau = 3;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'Full';
optl.Graph.nsz = [7,7];
optl.Graph.k = [];


imsz = optl.wsz+optl.psz-[1,1] ;
psz = optl.psz;
stpsz = [1,1];


Dversion = 'lenapatch+noise';
Wversion = 0;
eigop = 1;
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
    s = s+randn(size(s))*.1;
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
    s = s+randn(size(s))*.1;
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
    s = s+randn(size(s))*.1;
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
            if s_ref(i+psz(1)-1,j+psz(2)-1) == 1, % actually a line
                for k1 = -2:2
                    for k2 = -2:2
                        if i+k1>0 && j+k2>0 && i+k1<=optl.wsz(1) && j+k2<= optl.wsz(2)
                            if s_ref(i+k1+psz(1)-1,j+k2+psz(2)-1) == 1,
                                flind2 = rec2flat([i+k1,j+k2],imsz,psz,stpsz);
                                W(flind,flind2) = 1;
                            end
                        end
                    end
                end
            else
                for k1 = -1:1
                    for k2 = -1:1
                        if i+k1>0 && j+k2>0 && i+k1<=optl.wsz(1) && j+k2<= optl.wsz(2)
                            if s_ref(i+k1+psz(2)-1,j+k2+psz(2)-1) == 0,
                                flind2 = rec2flat([i+k1,j+k2],imsz,psz,stpsz);
                                W(flind,flind2) = 1;
                            end
                        end
                    end
                end
            end
        end
    end
    W = W+W';
    if strcmp(optl.Laplacian,'n')
        L = ulap(W);
    else
        L = nlap(W);
    end
    ind1 = [1,1];
    ind2 = optl.wsz;
    temp = L;
    L = {}; L{1}.ind1 = ind1;
    L{1}.ind2 = ind2; L{1}.M = temp;
    
end

% Version 0 Choose a scheme specified by optl
if Wversion == 0
    [L,scrop] = laplacian_from_image(s,optl);
    if(eigop)
        [phi,E] = eigs(L{1}.M,optl.neig,'sr');
        E = diag(E);
    end
end

%%%%%%%%%%%%%%%%%%%%% Visualizing and Testing Laplacian%%%%%%%%%%%%%%%%%%%%
% 
% 
%visualizing the eigenvector, eigenvalues
if(eigop)
    neig = optl.neig;
    Xeig = zeros(optl.wsz(1),optl.wsz(2),neig);
    for i = 1:neig
        foo = reshape(phi(:,i),optl.wsz(1),optl.wsz(2));
        Xeig(:,:,i) = foo';
    end
    o.colorbar = 0;
    square_plot(Xeig,o);
    figure;
    plot(E,'r');
end


figure;
a = L{1}.M;
for i = 1:size(a,1)
    a(i,i) = 0;
end
imagesc(-a);
colorbar;
title('Weight Matrix');

figure;
spy(a);


figure;
aa = sum(abs(a),2);
aa = reshape(aa,optl.wsz(1),optl.wsz(2))';
imagesc(aa);
colorbar;
title('Node Degree');

% 


% %%%%%%%%%%%%%%%%%%  Testing by Denoising Experiment %%%%%%%%%%%%%%%%%%%%%%%
% 
% % cbpdn with graph regularization
% mu = 4;
% lambda = .2;
% opt = {};
% opt.Verbose = 1;
% opt.MaxMainIter = 40;
% opt.rho = 100*lambda + 1;
% opt.RelStopTol = 2e-3;
% opt.AuxVarObj = 0;
% opt.HighMemSolve = 1;
% opt.L1Weight = 1;
% 
% 
% 
% % Load dictionary
% load([sporco_path '/Data/ConvDict.mat']);
% dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
% D = dmap('12x12x36');
% D = double(D);
% 
% opt.Ysolver = 'fista';
% [Xnl,~]= cbpdn_L(D,s,L,lambda,mu,opt);
% 
% [Xcn,~] = cbpdn(D,s,lambda,opt);
% 
% % reconstructing conv
% scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
%                                fft2(x)),3), 'symmetric');
% shrecnl = scnv(D,Xnl);
% figure; imagesc(sn); colorbar; title('noisy image'); colormap(gray);
% figure; imagesc(s); colorbar; title('noisy image high'); colormap(gray);
% 
% 
% figure; imagesc(shrecnl); colorbar; title('nonlocal high rec'); colormap(gray);
% figure; imagesc(shrecnl+sl); colorbar; title('nonlocal rec'); colormap(gray);
% figure; imagesc(sum3(abs(Xnl))); colorbar; title('nonlocal coef');colormap(gray);
% 
% 
% shrec = scnv(D,Xcn);
% 
% figure; imagesc(shrec); colorbar; title('conv high rec'); colormap(gray);
% figure; imagesc(shrec+sl); colorbar; title('conv rec'); colormap(gray);
% figure; imagesc(sum3(abs(Xcn))); colorbar; title('conv coef');colormap(gray);






