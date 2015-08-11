function Ind = patchind(iter,nsz,coefsz)
nc = coefsz(2);
k = 1;

rec_x = floor((iter-1)/nc)+1;
rec_y = iter-(rec_x-1)*nc;

l = min(nsz(1),rec_x-1);
r = min(nsz(1),coefsz(1) - rec_x);

u = min(nsz(2),rec_y-1);
d = min(nsz(2),nc-rec_y);

Ind = zeros(1,(r+l+1)*(d+u+1));

for i = -l:r
    for j = -u:d
        Ind(k) = iter+i*nc+j;
        k = k+1;
    end
end

