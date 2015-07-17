% split the photos into a larger training set

load('FlickrCC_512_512.mat');

nx = size(S,1);
ny = size(S,2);
m = size(S,3);
temp = zeros(nx/2,ny/2,m*4);


for i = 1:m
    temp(:,:,4*(i-1)+1) = S(1:nx/2,1:ny/2,i);
    temp(:,:,4*(i-1)+2) = S(1:nx/2,ny/2+1:ny,i);
    temp(:,:,4*(i-1)+3) = S(nx/2+1:nx,1:ny/2,i);
    temp(:,:,4*(i-1)+4) = S(nx/2+1:nx,ny/2+1:ny,i);
end
    

S = temp;
save('Flicker1_512_split','S');
