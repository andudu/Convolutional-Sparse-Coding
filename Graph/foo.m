% Test out distribution of eigenvalues


s = single(imresize(stdimage('lena.grey'),.5)) / 255;
[sl,sh] = lowpass(s,6,15);


wsz = [50,50];
psz = [8,8];
neig = 10;
tau = 1;

[L,scrop] = graphgen(sh, wsz, psz,neig,tau);

for i = 1:length(L)
    plot(L{i}.E);
    hold on;
end
hold off;


