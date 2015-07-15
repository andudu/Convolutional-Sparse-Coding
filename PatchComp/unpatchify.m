function A = unpatchify(X, patchSize, stepSize, sizeA)
% X is the flattened out patches, we want to reconstruct a matrix
% patchSize is the size of patches we took
% stepSize is the number of steps we take right and down for each patch
% sizeA is the size of the image you are trying to reproduce


m=sizeA(1);
n=sizeA(2);
k=1;
if size(sizeA,2)==3
    k=sizeA(3);
end

A=zeros(m,n,k);
Anormalizer=zeros(m,n,k);
i=1; j=1; counter=1;



while i+patchSize-1 <= m
    j=1;
    while j+patchSize-1<=n
        A(i:i+patchSize-1,j:j+patchSize-1,:)=A(i:i+patchSize-1,j:j+patchSize-1,:)+reshape(X(:,counter),patchSize,patchSize,k);
        Anormalizer(i:i+patchSize-1,j:j+patchSize-1,:)=Anormalizer(i:i+patchSize-1,j:j+patchSize-1,:)+ones(patchSize,patchSize,k);
        j=j+stepSize;
        counter=counter+1;
    end
    i=i+stepSize;
end

A=A./(Anormalizer+1e-4);

end