function [ ind ] = pnbr( xind, nsz,coefsz )
% Extract the rectangular index around a pixel neighborhood
% xind  = (i,j)  position of center 
% nsz = size of neighborhood (in each direction, e.g. nsz = 2 returns 9 points)
% imsz = xsz, ysz. use for boundary check. 

nc_max = coefsz(2);
nr_max = coefsz(1);

l1 = max(1,xind(2)-nsz); %left most patch
l2 = min(nc_max,xind(2)+nsz);

r1 = max(1, xind(1)-nsz);
r2 = min(nr_max,xind(1)+nsz);

a = l1:l2;
b = r1:r2;

a = a(ones(1,r2-r1+1),:);
a = reshape(a,1,numel(a));
b = repmat(b,1,l2-l1+1);
ind = [b;a];

end

