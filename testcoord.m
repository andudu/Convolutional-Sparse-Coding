s = single(rgbtogrey(stdimage('lena')))/255;
s = imresize(s,.5);

psz = [6,8];
st = [2,3];
imsz = size(s);
s1d = imageblocks(s,psz,st);
srec = zeros(size(s));
for i = 1:size(s1d,3)
    ind = flat2rec(i,imsz,psz,st);
    foo = rec2flat(ind,imsz,psz,st);
    srec(ind(1):ind(1)+psz(1)-1, ind(2):ind(2)+psz(2)-1) =   s1d(:,:,i);
end

figure;
imagesc(srec);