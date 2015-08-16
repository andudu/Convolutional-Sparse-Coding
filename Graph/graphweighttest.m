%scrit for testing 

%%%%%%%%%%%%%%%%%% Setting the Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setting options for Graph Laplacian
optl = {};
optl.wsz = [60,60];
optl.psz = [8,8];
optl.neig = 300;
optl.Lformat = 'Full';
optl.Laplacian = 'n';
optl.Graph.tau = 1;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'Window';
optl.Graph.nsz = [10,10];
optl.Graph.k = [];


imsz = optl.wsz+optl.psz-[1,1] ;
psz = optl.psz;
stpsz = [1,1];
load('stdnoise.mat');
n = r_noise(1:imsz(1),1:imsz(2));
Dversion = 'simpleline+noise';
Wversion = 1;
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
% Perfect graph for the simple line case
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
    
end

% % Version 0 Choose a scheme specified by optl
% if Wversion == 0
%     [L,scrop] = laplacian_from_image(s,optl);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% The Actual Test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% a = L{1}.M;
% for i = 1:size(a,1)
%     a(i,i) = 0;
% end
% 
% %test1;
% optw.imsz = imsz ;
% optw.psz = psz;
% optw.BackgroundImage = s;
% optw.W = -a;
% ind1 = [35,45]; % a curved line in lena patch
% ind1 = [20, 20];   %a point on the simple line (long line)
% [s1,~] = showweights(ind1,optw);
% 
% 
% optw = {};
% optw.imsz = imsz ;
% optw.psz = psz;
% optw.BackgroundImage = s;
% optw.tau = optl.Graph.tau;
% optw.Metric = optl.Graph.Metric;
% optw.GraphType = optl.Graph.GraphType;
% optw.nsz = optl.Graph.nsz;
% optw.coefsz = imsz - psz + [1,1];
% [s2,~] = showweights(ind1,optw);
% 
% 
% figure;
% subplot(1,2,1);
% imagesc(s1);
% hold on;
% plot(ind1(2),ind1(1),'rx');
% hold off;
% title('Normalized Weights');
% 
% subplot(1,2,2);
% imagesc(s2);
% hold on;
% plot(ind1(2),ind1(1),'rx');
% hold off;
% title('Original Weights');
% 
% 
% 
% ind2 = [24,22];  % a "confetti patch"
% optw.W = -a;
% [s1,~] = showweights(ind2,optw);
% optw = rmfield(optw,'W');
% [s2,~] = showweights(ind2,optw);
% figure;
% subplot(1,2,1);
% imagesc(s1);
% hold on;
% plot(ind2(2),ind2(1),'rx');
% hold off;
% title('Normalized Weights');
% 
% subplot(1,2,2);
% imagesc(s2);
% hold on;
% plot(ind2(2),ind2(1),'rx');
% hold off;
% title('Original Weights');


%test3;  % full image
optw = {};
s = double(stdimage('lena.grey')) / 255;
s = imresize(s,.5);
s = s+r_noise;
[sl,sh] = lowpass(s,6,15);
imsz = size(sh);
psz = [12,12];
s = sh;
ind1 = [53,106]; %This lies on one of the edges
%ind2 = [68,213]; % This looks like a confetti patch
optw.imsz = imsz ;
optw.psz = psz;
optw.BackgroundImage = s;
optw.tau = 2;
optw.Metric = 'Cosine';
optw.GraphType = 'Window';
optw.nsz = [15,15];
optw.coefsz = imsz - psz + [1,1];

[s1,~] = showweights(ind1,optw);
figure; 
imagesc(s1);
hold on;
plot(ind1(2),ind1(1),'rx');
hold off;
title('Original Weights');

% [s2,~] = showweights(ind2,optw);
% figure; 
% imagesc(s2);
% hold on;
% plot(ind2(2),ind2(1),'rx');
% hold off;
% title('Original Weights');



