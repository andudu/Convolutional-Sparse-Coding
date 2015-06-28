function Y = convsum(X,D,ind)
%find the partial convolutional sum correspondint to ind
Y = zeros(size(X,1),size(X,2));
for i = ind
    Y = Y + ifft2(bsxfun(@times,fft2(D(:,:,i),size(X,1),size(X,2)),...
        fft2(X(:,:,i))));
end
