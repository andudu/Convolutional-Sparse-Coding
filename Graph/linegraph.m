%script for testing simple graphs on an image of a line.

%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imsz = [69,69];
psz = [12,12];
coefsz = imsz-psz+[1,1];
stpsz = [1,1];

Wversion = 2;
Lversion = 'n';
Dversion = 'simpleline';

%%%%%%%%%%%%%%%%%%   Generating Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Dversion, 'simpleline'),
    % simple line
    s = zeros(69,69);
    s(35:36,:) = 1;
    s_ref = s;
    indfoo = 34*60+20;
end

if strcmp(Dversion, 'simpleline+noise'),
    % simple line
    s = zeros(69,69);
    s(35:36,:) = 1;
    s = s+randn(size(s))*.1;
    s_ref = s;
    indfoo = 34*60+20;
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
end

% Version 2. Do full nonlocal graph
if Wversion == 2,    
    tau = 15;
    s1d = imageblocks(s,psz,stpsz); 
    s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));
    W = sqdist(s1d,s1d);
    W = exp(-W/(tau));
end

% Version 3. K nearest Neighbors of full graph
if Wversion == 3, 
    tau = 1;
    k = 40;
    s1d = imageblocks(s,psz,stpsz); 
    s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));
    W = sqdist(s1d,s1d);
    W = exp(-W/(tau));
    for iter = 1:size(W,2)
        [~,is] = sort(W(:,iter),'descend');
        W(is(k+1:end),iter) = 0;
    end
    W = max(W,W');   
end


% Version 4. K nearest Neighbors within a smaller window
if Wversion == 4, 
    tau = 1;
    k = 10;
    nsz = [3,3];
    s1d = imageblocks(s,psz,stpsz); 
    s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));
    W = sqdist(s1d,s1d);
    W = exp(-W/(tau));
    for iter = 1:size(W,1);
        Ind = patchind(iter,nsz,coefsz);
        [v,is] = sort(W(Ind,iter),'descend');
        W(:,iter) = 0;
        W(Ind(is(1:k)),iter) = v(1:k);
    end
    W = max(W,W');          
end







%%%%%%%%%%%%%%%%%%%%% Visualizing and Testing Laplacian%%%%%%%%%%%%%%%%%%%%

if Lversion == 'u'
    L = ulap(W);
end
if Lversion == 'n'
    L = nlap(W);
end

%visualizing the eigenvector, eigenvalues

neig = 40;
[phi,E] = eigs(L,neig,'sr'); E = diag(E);

Xeig = zeros(coefsz(1),coefsz(2),neig);
for i = 1:neig
    foo = reshape(phi(:,i),coefsz(1),coefsz(2));
    Xeig(:,:,i) = foo';
end
o.colorbar = 1;
square_plot(Xeig,o);

figure;
plot(E,'r');

figure;
imagesc(W);
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
 





