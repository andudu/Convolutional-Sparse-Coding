function M = spmask(imsz,nsz)

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
return;

