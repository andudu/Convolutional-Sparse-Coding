% scrip for generating eigenvectors / sparse graphs

%%%%%%%%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standard Training images
S0 = zeros(512, 512, 5, 'single');
S0(:,:,1) = single(stdimage('lena.grey')) / 255;
S0(:,:,2) = single(stdimage('barbara.grey')) / 255;
S0(:,:,3) = single(stdimage('kiel.grey')) / 255;
S0(:,:,4) = single(rgb2gray(stdimage('mandrill'))) / 255;
tmp = single(stdimage('man.grey')) / 255;
S0(:,:,5) = tmp(101:612, 101:612);


%Reduce images size to speed up demo script
% tmp = [];
% for k = 1:size(S0,3),
%   tmp(:,:,k) = imresize(S0(:,:,k), 0.5);
% end
% S0 = tmp;

tag1 = 'Standard512';

% % Flicker Training Images
% tag1 = 'Flicker';
% load('Flicker1_512_split.mat');
% S0 = S;
% S0 = double(S0)/255;
% clear S;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Eigenvectors %%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(['CacheData/',tag1],'dir')
    mkdir(['CacheData/',tag1]);
end


imnum = size(S0,3);

lambda = 4;


% % Nystrom (currently not working for some reason...)
% for i = 1:imnum
%     i
%     temp = double(S0(:,:,i))/255;
%     [~,sh] = lowpass(temp,lambda,15);
%     clear temp;
%     imsz = size(sh);
%     psz = [12,12];
%     stpsz = [1,1];
%     sh = padarray(sh,psz-[1,1],'symmetric','post');
%     
%     %generate graph
%     opt.tau = 1.2;
%     opt.Laplacian = 'n';
%     opt.numsample = 800;
%     opt.Metric = 'Euclidean';
%     opt.neig = 30;
%     data = imageblocks(sh,psz,stpsz);
%     data = reshape(data,size(data,1)*size(data,2),size(data,3));
%     data = data';
%     [phi,E] = nystrom(data,opt);
%     save(['CacheData/',foldname, '/eig',num2str(i),'.mat'],'phi','E');
%     
% end

% The Eigenvectors

for i = 5:imnum
    i
    temp = S0(:,:,i);
    [~,sh] = lowpass(temp,lambda,15);
    clear temp;
    imsz = size(sh);
    psz = [12,12];
    stpsz = [1,1];
    sh = padarray(sh,psz-[1,1],'symmetric','post');
    
    %generate graph
    optl = {};
    optl.wsz = imsz;
    optl.psz = [12,12];
    optl.neig = 100;
    optl.Lformat = 'Sparse';
    optl.Graph.Laplacian = 'n';
    optl.Graph.tau = 1;
    optl.Graph.Metric = 'Euclidean';
    optl.Graph.GraphType = 'KNearest';
    optl.SaveMem = 1;
    % optl.Threshold = 1e-4;
    optl.Graph.nsz = [];
    optl.Graph.k = 30;
    disp('generating graph');
    tic;
    [L,~] = laplacian_from_image(sh,optl);
    toc
    disp('calculating eig');
    tic;
    [phi,E] = eigs(L{1}.M,optl.neig,'sr');
    toc
    E = diag(E);
    
    %         temp = reshape(phi,256,256,optl.neig);
    %         temp = permute(temp,[2,1,3]);
    %         square_plot(temp(:,:,50:60),{});
    %         figure; plot(E);
    
    if ~exist(['CacheData/',tag1,'/Eig'],'dir')
        mkdir(['CacheData/',tag1,'/Eig']);
    end
    phi = single(phi);
    save(['CacheData/',tag1, '/Eig/',optl.Graph.Metric,'eig',num2str(i),'.mat'],'phi','E');
    
    if ~exist(['CacheData/',tag1,'/Mat'],'dir')
        mkdir(['CacheData/',tag1,'/Mat']);
    end
    save(['CacheData/',tag1,'/Mat', '/',optl.Graph.Metric,'M',num2str(i),'.mat'],'L');    
    clear L phi E;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imsz = [size(S0,1),size(S0,2)];
save(['CacheData/',tag1,'/',optl.Graph.Metric,'param.mat'],'imnum','lambda','imsz','psz','optl');




% %%% Next Batch %%%
% % scrip for generating eigenvectors / sparse graphs
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Standard Training images
% S0 = zeros(512, 512, 5, 'single');
% S0(:,:,1) = single(stdimage('lena.grey')) / 255;
% S0(:,:,2) = single(stdimage('barbara.grey')) / 255;
% S0(:,:,3) = single(stdimage('kiel.grey')) / 255;
% S0(:,:,4) = single(rgb2gray(stdimage('mandrill'))) / 255;
% tmp = single(stdimage('man.grey')) / 255;
% S0(:,:,5) = tmp(101:612, 101:612);
% 
% 
% % %Reduce images size to speed up demo script
% % tmp = [];
% % for k = 1:size(S0,3),
% %   tmp(:,:,k) = imresize(S0(:,:,k), 0.5);
% % end
% % S0 = tmp;
% 
% % % Flicker Training Images
% % tag1 = 'Flicker';
% % load('Flicker1_512_split.mat');
% % S0 = S;
% % S0 = double(S0)/255;
% % clear S;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Eigenvectors %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if ~exist(['CacheData/',tag1],'dir')
%     mkdir(['CacheData/',tag1]);
% end
% 
% 
% imnum = size(S0,3);
% 
% lambda = 4;
% 
% 
% % % Nystrom (currently not working for some reason...)
% % for i = 1:imnum
% %     i
% %     temp = double(S0(:,:,i))/255;
% %     [~,sh] = lowpass(temp,lambda,15);
% %     clear temp;
% %     imsz = size(sh);
% %     psz = [12,12];
% %     stpsz = [1,1];
% %     sh = padarray(sh,psz-[1,1],'symmetric','post');
% %     
% %     %generate graph
% %     opt.tau = 1.2;
% %     opt.Laplacian = 'n';
% %     opt.numsample = 800;
% %     opt.Metric = 'Euclidean';
% %     opt.neig = 30;
% %     data = imageblocks(sh,psz,stpsz);
% %     data = reshape(data,size(data,1)*size(data,2),size(data,3));
% %     data = data';
% %     [phi,E] = nystrom(data,opt);
% %     save(['CacheData/',foldname, '/eig',num2str(i),'.mat'],'phi','E');
% %     
% % end
% 
% % The Eigenvectors
% 
% for i = 1:imnum
%     i
%     temp = S0(:,:,i);
%     [~,sh] = lowpass(temp,lambda,15);
%     clear temp;
%     imsz = size(sh);
%     psz = [12,12];
%     stpsz = [1,1];
%     sh = padarray(sh,psz-[1,1],'symmetric','post');
%     
%     %generate graph
%     optl = {};
%     optl.wsz = imsz;
%     optl.psz = [12,12];
%     optl.neig = 100;
%     optl.Lformat = 'Sparse';
%     optl.Graph.Laplacian = 'n';
%     optl.Graph.tau = 1;
%     optl.Graph.Metric = 'Cosine';
%     optl.Graph.GraphType = 'KNearest';
%     optl.SaveMem = 1;
%     % optl.Threshold = 1e-4;
%     optl.Graph.nsz = [];
%     optl.Graph.k = 30;
%     disp('generating graph');
%     tic;
%     [L,~] = laplacian_from_image(sh,optl);
%     toc
%     disp('calculating eig');
%     tic;
%     [phi,E] = eigs(L{1}.M,optl.neig,'sr');
%     toc
%     E = diag(E);
%     
%     %         temp = reshape(phi,256,256,optl.neig);
%     %         temp = permute(temp,[2,1,3]);
%     %         square_plot(temp(:,:,50:60),{});
%     %         figure; plot(E);
%     
%     if ~exist(['CacheData/',tag1,'/Eig'],'dir')
%         mkdir(['CacheData/',tag1,'/Eig']);
%     end
%     phi = single(phi);
%     save(['CacheData/',tag1, '/Eig/',optl.Graph.Metric,'eig',num2str(i),'.mat'],'phi','E');
%     
%     if ~exist(['CacheData/',tag1,'/Mat'],'dir')
%         mkdir(['CacheData/',tag1,'/Mat']);
%     end
%     save(['CacheData/',tag1,'/Mat', '/',optl.Graph.Metric,'M',num2str(i),'.mat'],'L');    
%     clear L phi E;
%     
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imsz = [size(S0,1),size(S0,2)];
% save(['CacheData/',tag1,'/',optl.Graph.Metric,'param.mat'],'imnum','lambda','imsz','psz','optl');






