function y = PxT(A)
y = zeros(2*size(A,1),2*size(A,2),size(A,3),size(A,4));
for i = 1:size(A,3)
    for j = 1:size(A,4)
        y(:,:,i,j) = kron(A(:,:,i,j),[1,0;0,0]);
    end
end
return