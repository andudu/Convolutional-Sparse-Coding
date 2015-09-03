% script for testing the masking matrix 

imsz = 12;
nsz = 2;
ind = 34;

I = 1:1:imsz;
J = 1:1:imsz;

for k = 1:nsz-1
    I = [I,1:1:imsz-k];
    J = [J,k+1:1:imsz];
end

for k = 1:nsz-1
    I = [I,k+1:1:imsz];
    J = [J,1:1:imsz-k];
end

r = ones(size(I));

M = sparse(I,J,r);

M = kron(M,M);

for ind = 1:1:imsz*imsz

x = zeros(imsz*imsz,1);
x(ind,1) = 1;

y = M*x;
imagesc(reshape(y,imsz,imsz));
drawnow;
pause(.3);

end





    


    