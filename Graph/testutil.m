%script for testing laplacian_from_image

s = single(stdimage('lena.grey')) / 255;
s = imresize(s,.5);

%setting options for Graph Laplacian
optl.wsz = [60,60];
optl.psz = [12,12];
optl.neig = 30;
optl.Lformat = 'Eig';
optl.Laplacian = 'n';
optl.Graph.tau = 3;
optl.Graph.Metric = 'Euclidean';
optl.Graph.GraphType = 'WindowKNearest';
optl.Graph.nsz = [4,4];
optl.Graph.k = 5;


%generating graph laplacian
[sl,sh] = lowpass(s,5,15);
[L,shcrop] = laplacian_from_image(sh,optl);


%plotting eigenvalues and vectors
X = [];
for i = 1:length(L)
   X(:,:,i) = reshape(L{i}.phi(:,2),optl.wsz(1),optl.wsz(2))'; 
end
square_plot(X,{});
figure;
for i = 1:length(L)
    plot(L{i}.E);
    hold on;
end
hold off;


% 
% %plotting the weight matrix
% X = [];
% for i = 1:length(L)
%    X(:,:,i) = L{i}.M; 
% end
% square_plot(X,{});
