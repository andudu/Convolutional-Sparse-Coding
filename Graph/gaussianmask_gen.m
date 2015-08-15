function W = gaussianmask_gen(nsz, coefsz)
% generate a gaussian mask for every point in coefsz 
% and build a matrix

gsz = nsz*2+[1,1];
m = gauss2d(gsz,gsz(1)/4);
W = zeros(coefsz(1)*coefsz(2),coefsz(1)*coefsz(2));
for i = 1:coefsz(1);
    for j = 1:coefsz(2);
        a = max(1,i-nsz(1));
        b = max(1,j-nsz(2));
        c = min(coefsz(1), i+nsz(1));
        d = min(coefsz(1), j + nsz(2));
        shifti = i-(nsz(1)+1);
        shiftj = j-(nsz(2)+1);
        mtemp = m(a-shifti:c-shifti,b-shiftj:d-shiftj);
        k = (i-1)*(coefsz(2))+j;
        shiftk = (i-a)*(d-b+1)+j-1;
        mtemp = reshape(mtemp',size(mtemp,1)*size(mtemp,2),1);
        W(k,k-shiftk:k-shiftk+length(mtemp)-1) = mtemp;         
    end
    W = max(W,W');
end

W = sparse(W);

return;