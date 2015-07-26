%testing nonlocal lowpass

%load precomputed graph info

load('data_graph.mat');
s = sl+sh;

foo = sh+.1*randn(size(sh));
[a,b] = lowpass(foo,10);


for i = 1:2:21
    [sln,slh] = lowpass_nl(foo,L,i);
    figure;
    imagesc(sln);
    colorbar;
end



figure;
imagesc(a);
colorbar;



for i = 1:length(L)
    plot(L{i}.E);
    hold on;
end
hold off;