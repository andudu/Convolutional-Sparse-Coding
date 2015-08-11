%script for testing simple graphs on an image of a line.

%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting options for Graph Laplacian
optl = {};
optl.wsz = [60,60];
optl.psz = [12,12];
optl.neig = 30;
optl.Lformat = 'Full';
optl.Laplacian = 'n';
optl.Graph.tau = .8;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'Full';
optl.Graph.nsz = [5,5];
optl.Graph.k = [];


imsz = optl.wsz+optl.psz-[-1,-1] ;
psz = optl.psz;
coefsz = imsz-psz+[1,1];
stpsz = [1,1];


Dversion = 'lenapatch+noise';
Wversion = 0;
%%%%%%%%%%%%%%%%%%   Generating Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Dversion, 'simpleline'),
    % simple line
    s = zeros(imsz);
    s(35:36,:) = 1;
    s_ref = s;
end

if strcmp(Dversion, 'simpleline+noise'),
    % simple line
    s = zeros(imsz);
    s(35:36,:) = 1;
    s_ref = s;
    s = s+randn(size(s))*.2;
end

if strcmp(Dversion, 'lenapatch'),
    %patch from lena
    s = single(stdimage('lena.grey')) / 255;
    s = imresize(s,.5);
    [sl,sh] = lowpass(s,7,15);
    s = sh(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
end

if strcmp(Dversion, 'lenapatch+noise'),
    %patch from lena
    s = single(stdimage('lena.grey')) / 255;
    s = imresize(s,.5);
    [sl,sh] = lowpass(s,7,15);
    s = sh(50:1:imsz(1)+50-1,160:1:imsz(2)+160-1);
    s = randn(size(s))*.1+s;
end


%%%%%%%%%%%%%%%%%   Building the Graph Weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version 1. Connect nonzero/zero elements to its nonzero/zero neighbors
% Perfect paring in this case
if Wversion == 1
    W = sparse(coefsz(1)*coefsz(2), coefsz(1)*coefsz(2));
    for i = 1:coefsz(1)
        for j = 1:coefsz(2)
            flind = rec2flat([i,j],imsz,psz,stpsz);
            if s_ref(i,j) == 1, % actually a line
                for k1 = -2:2
                    for k2 = -2:2
                        if i+k1>0 && j+k2>0 && i+k1<=coefsz(1) && j+k2<= coefsz(2)
                            if s_ref(i+k1,j+k2) == 1,
                                flind2 = rec2flat([i+k1,j+k2],imsz,psz,stpsz);
                                W(flind,flind2) = 1;
                            end
                        end
                    end
                end
            else
                for k1 = -1:1
                    for k2 = -1:1
                        if i+k1>0 && j+k2>0 && i+k1<=coefsz(1) && j+k2<= coefsz(2)
                            if s_ref(i+k1,j+k2) == 0,
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
end

% Version 0 Choose a scheme specified by optl
if Wversion == 0
    [L,scrop] = laplacian_from_image(s,optl);
    [phi,E] = eigs(L{1}.M,optl.neig,'sr');
    E = diag(E);
end

%%%%%%%%%%%%%%%%%%%%% Visualizing and Testing Laplacian%%%%%%%%%%%%%%%%%%%%


%visualizing the eigenvector, eigenvalues

neig = optl.neig;
Xeig = zeros(optl.wsz(1),optl.wsz(2),neig);
for i = 1:neig
    foo = reshape(phi(:,i),optl.wsz(1),optl.wsz(2));
    Xeig(:,:,i) = foo';
end
o.colorbar = 1;
square_plot(Xeig,o);

figure;
plot(E,'r');

figure;
a = L{1}.M;
for i = 1:size(a,1)
    a(i,i) = 0;
end
imagesc(-a);
colorbar;





% %%%%%%%%%%%%%%%%%%  Testing by Denoising Experiment %%%%%%%%%%%%%%%%%%%%%%%
% 
% % % cbpdn with graph regularization
% mu = 10;
% lambda = 0.2;
% opt = {};
% opt.Verbose = 1;
% opt.MaxMainIter = 130;
% opt.rho = 100*lambda + 1;
% opt.RelStopTol = 1e-3;
% opt.AuxVarObj = 0;
% opt.HighMemSolve = 1;
% opt.L1Weight = 1;
% temp = L;
% L= {};
% L{1}.M = temp;
% L{1}.ind1 = [1,1];
% L{1}.ind2 = coefsz;
% 
% % Load dictionary
% load([sporco_path '/Data/ConvDict.mat']);
% dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
% D = dmap('12x12x36');
% D = double(D);
% s = s(1:coefsz(1),1:coefsz(2));
% 
% [Xcn,~] = cbpdn(D,s,lambda,opt);
% disp('cbpdn opt');
% 
% opt.Ysolver = 'fista';
% lambda = .1;
% [Xnl,~]= cbpdn_L(D,s,L,lambda,mu,opt);
% disp('nl opt');
% 
% % reconstructing conv
% scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
%                                fft2(x)),3), 'symmetric');
% shrecnl = scnv(D,Xnl);
% shreccn = scnv(D,Xcn);
% 
% figure; imagesc(shrecnl); colorbar; title('nonlocal rec'); colormap(gray);
% figure; imagesc(sum3(abs(Xnl))); colorbar; title('nonlocal coef');colormap(gray);
% 
% figure; imagesc(shreccn); colorbar; title('conv rec');colormap(gray);
% figure; imagesc(sum3(abs(Xcn))); colorbar; title('conv coef');colormap(gray);
 





