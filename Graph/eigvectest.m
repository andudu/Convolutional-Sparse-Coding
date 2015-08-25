% eigvector test: 

% Testing effects of Graph Parameters on Eigenvectors and Eigenvalues


%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting options for Graph Laplacian
optl = {};
optl.wsz = [60,60];
optl.psz = [8,8];
optl.neig = 30;
optl.Lformat = 'Full';
optl.Laplacian = 'n';
optl.Graph.tau = 1;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'Full';
optl.Graph.nsz = [8,8];
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

% Version 0 Choose a scheme specified by optl

optl.Graph.nsz = [8,8];

if Wversion == 0
    [L1,scrop] = laplacian_from_image(s,optl);
    if(eigop)
        [phi1,E1] = eigs(L1{1}.M,optl.neig,'sr');
        E1 = diag(E1);
    end
end



%%%%%%%%%%%%%%%%%%%%% Visualizing and Testing Laplacian%%%%%%%%%%%%%%%%%%%%
% 
% Plot 3 eigenvectors


%visualizing the eigenvector, eigenvalues
Xeig1 = [];
for i = 2:10
    foo = reshape(phi1(:,i),optl.wsz(1),optl.wsz(2));
    Xeig1(:,:,i-1) = foo';
end
o.colorbar = 0;
square_plot(Xeig1,o);
figure; 
plot(E1);


%generate graph
opt = {};
opt.tau = 1;
opt.Laplacian = 'n';
opt.numsample = 500;
opt.Metric = 'Euclidean';
opt.neig = 10;


data = imageblocks(s,psz,stpsz);
data = reshape(data,size(data,1)*size(data,2),size(data,3));
data = data';
[phi,E] = nystrom(data,opt);
disp('Normalized graph generated');
X1 = [];
for i = 2:opt.neig,
    v = phi(:,i);
    v = reshape(v,[60,60])';
    X1(:,:,i-1) = v;
end
square_plot(X1,{});
figure; 
plot(E);



% a = L{1}.M;
% for i = 1:size(a,1)
%     a(i,i) = 0;
% end
% figure;
% imagesc(-a);
% colorbar;
% title('Weight Matrix');
% 
% figure;
% spy(a);
% 
% figure;
% aa = sum(abs(a),2);
% aa = reshape(aa,optl.wsz(1),optl.wsz(2))';
% imagesc(aa);
% colorbar;
% title('Node Degree');
