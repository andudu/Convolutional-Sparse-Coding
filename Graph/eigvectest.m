% eigvector test: 
% Testing effects of Graph Parameters on Eigenvectors and Eigenvalues

%script for testing simple graphs on an image of a line.

%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting options for Graph Laplacian
optl = {};
optl.wsz = [60,60];
optl.psz = [8,8];
optl.neig = 300;
optl.Lformat = 'Sparse';
optl.Laplacian = 'n';
optl.Graph.tau = 2;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'Window';
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

optl.Graph.nsz = [3,3];

if Wversion == 0
    [L1,scrop] = laplacian_from_image(s,optl);
    if(eigop)
        [phi1,E1] = eigs(L1{1}.M,optl.neig,'sr');
        E1 = diag(E1);
    end
end

optl.Graph.nsz = [5,5];

% Version 0 Choose a scheme specified by optl
if Wversion == 0
    [L2,scrop] = laplacian_from_image(s,optl);
    if(eigop)
        [phi2,E2] = eigs(L2{1}.M,optl.neig,'sr');
        E2 = diag(E2);
    end
end

optl.Graph.nsz = [9,9];

% Version 0 Choose a scheme specified by optl
if Wversion == 0
    [L3,scrop] = laplacian_from_image(s,optl);
    if(eigop)
        [phi3,E3] = eigs(L3{1}.M,optl.neig,'sr');
        E3 = diag(E3);
    end
end


%%%%%%%%%%%%%%%%%%%%% Visualizing and Testing Laplacian%%%%%%%%%%%%%%%%%%%%
% 
% Plot 3 eigenvectors
x = 1:1:optl.neig;
plot(x,E1,x,E2,x,E3);
legend('nsz = 3', 'nsz = 5' , 'nsz = 9');
saveas(gcf,'nsz-eigvalue','png');

%visualizing the eigenvector, eigenvalues
Xeig1 = [];
for i = 15:29
    foo = reshape(phi1(:,i),optl.wsz(1),optl.wsz(2));
    Xeig1(:,:,i-14) = foo';
end
o.colorbar = 0;
square_plot(Xeig1,o);
saveas(gcf,'nsz3-eigvector','png');

Xeig2 = [];
for i = 15:29
    foo = reshape(phi2(:,i),optl.wsz(1),optl.wsz(2));
    Xeig2(:,:,i-14) = foo';
end
o.colorbar = 0;
square_plot(Xeig2,o);
saveas(gcf,'nsz5-eigvector','png');

Xeig3 = [];
for i = 15:29
    foo = reshape(phi3(:,i),optl.wsz(1),optl.wsz(2));
    Xeig3(:,:,i-14) = foo';
end
o.colorbar = 0;
square_plot(Xeig3,o);
saveas(gcf,'nsz9-eigvector','png');


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
