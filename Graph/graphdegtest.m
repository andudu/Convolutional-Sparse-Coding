%graph connectivity test

%script for testing and plotting cbpdn_L


%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load  Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load a saved noise
imsz = [256,256];
snoise = randn(256,256)*.1;
sref = double(stdimage('lena.grey')) / 255;
sref = imresize(sref,.5);
s = sref+snoise;
[sl,sh] = lowpass(s,5,15);


%generate graph
optl = {};
optl.wsz = [60,60];
optl.psz = [12,12];
optl.neig = 60;
optl.Lformat = 'Sparse';
optl.Laplacian = 'u';
optl.Graph.tau = 1;
optl.Graph.Metric = 'Cosine';
optl.Graph.GraphType = 'Window';
optl.Graph.nsz = [9,9];
optl.Graph.k = [];
[L,sh] = laplacian_from_image(sh,optl);
disp('graph generated');



%cropping the image
sl = sl(1:size(sh,1),1:size(sh,2));
sref = sref(1:size(sh,1),1:size(sh,2));
[~,shref] = lowpass(sref,5,15);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xdeg = zeros(size(sh));
wsz = optl.wsz;
nwpr = floor(imsz(2)/wsz(2));
for k = 1:length(L)
    temp = L{k}.M;
    temp = full(temp);
    num = zeros(1,size(temp,1));
    for j = 1:size(temp,1)
        temp(j,j) = 0;
        foo = temp(j,:);
        num(j) = numel(foo(foo>1e-4));
        if num(j) == 0
            num(j) = 1;
        end
    end
    a = sum(-temp,1);
    a = a./num;
    a = reshape(a,wsz(1),wsz(2))';
    indi = floor((k-1)/nwpr)+1;
    indj = k - nwpr*(indi-1);
    indx = (indi-1)*wsz(1) + 1;
    indy = (indj-1)*wsz(2) + 1;
    xdeg(indx:indx + wsz(2)-1, indy:indy + wsz(1)-1) = a;
    clear temp;
end

figure; imagesc(xdeg); colorbar;



