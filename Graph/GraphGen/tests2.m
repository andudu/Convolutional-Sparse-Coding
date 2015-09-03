% test for testing the WindowKNearest functionality of gen


% % Standard Training images
S0 = double(stdimage('lena.grey')) / 255;
S0 = imresize(S0,.25);


% Building the Graph


[~,sh] = lowpass(S0,5,15);
imsz = size(sh);
psz = [12,12];
stpsz = [1,1];

%generate graph
optl = {};
optl.wsz = imsz;
optl.psz = [12,12];
optl.neig = [];
optl.Lformat = 'Sparse';
optl.Graph.Laplacian = 'n';
optl.Graph.tau = 1;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'WindowKNearest';
optl.SaveMem = 1;
optl.Graph.nsz = [5,5];
optl.Graph.k = 20;

disp('generating Cosine graph');
[L,~] = laplacian_from_image(sh,optl);





    

