%script for testing simple graphs on an image of a line.

%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imsz = [69,69];
psz = [10,10];
coefsz = imsz-psz+[1,1];
stpsz = [1,1];

Wversion = 3;
Lversion = 'n';
Dversion = 'simpleline+noise';

%%%%%%%%%%%%%%%%%%   Generating Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Dversion, 'simpleline'),
    % simple line
    s = zeros(69,69);
    s(35:36,:) = .5;
    indfoo = 34*60+20;
end

if strcmp(Dversion, 'simpleline+noise'),
    % simple line
    s = zeros(69,69);
    s(35:36,:) = .5;
    s = s+randn(size(s))*.1;
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
            if s(i,j) == 1, % actually a line
                for k1 = -2:2
                    for k2 = -2:2
                        if i+k1>0 && j+k2>0 && i+k1<=coefsz(1) && j+k2<= coefsz(2)
                            if s(i+k1,j+k2) == 1,
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
                            if s(i+k1,j+k2) == 0,
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
    np = psz(1)*psz(2);     
    tau = 1;
    s1d = imageblocks(s,psz,stpsz); 
    s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));
    W = sqdist(s1d,s1d);
    W = exp(-W/(np*tau));
end

% Version 3. K nearest Neighbors of full graph
if Wversion == 3,
    np = psz(1)*psz(2);    
    tau = 1;
    k = 10;
    s1d = imageblocks(s,psz,stpsz); 
    s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));
    W = sqdist(s1d,s1d);
    W = exp(-W/(tau*np));
    for iter = 1:size(W,1)
        [~,is] = sort(W(:,iter),'descend');
        W(is(k+1:end),iter) = 0;
    end
    W = max(W,W');
    W = sparse(W);
end


% Version 4. K nearest Neighbors within a smaller window
if Wversion == 4
    np = psz(1)*psz(2);
    tau = 1;
    k = 6;
    nsz = [2,2];
    s1d = imageblocks(s,psz,stpsz); 
    s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));
    W = sqdist(s1d,s1d);
    W = exp(-W/(tau*np));
    for iter = 1:size(W,1);
        Ind = patchind(iter,nsz,coefsz);
        [v,is] = sort(W(Ind,iter),'descend');
        W(:,iter) = 0;
        W(Ind(is(1:k)),iter) = v(1:k);
    end
    W = max(W,W');
    W = sparse(W);    
        
end

% Version 5. Adding a spatial weighting term to original matrix






%%%%%%%%%%%%%%%%%%%%% Visualizing and Testing Laplacian%%%%%%%%%%%%%%%%%%%%

if Lversion == 'u'
    L = ulap(W);
end
if Lversion == 'n'
    L = nlap(W);
end
nnz(L)/numel(L)

neig = 40;
[phi,E] = eigs(L,neig,'sr'); E = diag(E);

Xeig = zeros(coefsz(1),coefsz(2),neig);
for i = 1:neig
    foo = reshape(phi(:,i),coefsz(1),coefsz(2));
    Xeig(:,:,i) = foo';
end

square_plot(Xeig,{});

figure;
plot(E,'r');



