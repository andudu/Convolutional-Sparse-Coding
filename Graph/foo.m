%script for testing simple graphs on an image of a line.

%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting options for Graph Laplacian
optl = {};
optl.wsz = [60,60];
optl.psz = [8,8];
optl.neig = 300;
optl.Lformat = 'Full';
optl.Laplacian = 'u';
optl.Graph.tau = 3;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'Full';
optl.Graph.nsz = [8,8];
optl.Graph.k = [];


imsz = optl.wsz+optl.psz-[1,1] ;
psz = optl.psz;
stpsz = [1,1];


Dversion = 'lenapatch+noise';
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
end


a = L{1}.M;
for i = 1:size(a,1)
    a(i,i) = 0;
end

%test1;
optw.imsz = imsz ;
optw.psz = psz;
optw.BackgroundImage = s;
optw.W = -a;
ind = [35,45]; % a line
[s1,s2] = showweights(ind,optw);
figure; 
imagesc(s1);
title('Normalized Weights');

%test2;
optw = {};
[sl,sh] = lowpass(s,6,15);
imsz = size(sh);
psz = [12,12];
s = sh;
% ind = [51,158]; %This lies on one of the edges
ind = [35,45]; % This looks like a confetti patch
optw.imsz = imsz ;
optw.psz = psz;
optw.BackgroundImage = s;
optw.tau = 3;
optw.Metric = 'Cosine';
optw.GraphType = 'Full';
optw.coefsz = imsz - psz + [1,1];

[s1,s2] = showweights(ind,optw);
figure; 
imagesc(s1);
title('weights');



% %test2;
% optw = {};
% s = double(stdimage('lena.grey')) / 255;
% s = imresize(s,.5);
% s = s+.1*randn(size(s));
% [sl,sh] = lowpass(s,6,15);
% imsz = size(sh);
% psz = [12,12];
% s = sh;
% % ind = [51,158]; %This lies on one of the edges
% ind = [63,47]; % This looks like a confetti patch
% optw.imsz = imsz ;
% optw.psz = psz;
% optw.BackgroundImage = s;
% optw.tau = 2;
% optw.Metric = 'Cosine';
% optw.GraphType = 'Full';
% optw.nsz = [12,12];
% optw.coefsz = imsz - psz + [1,1];
% 
% [s1,s2] = showweights(ind,optw);
% figure; 
% imagesc(s1);
% figure;
% imagesc(s2);




