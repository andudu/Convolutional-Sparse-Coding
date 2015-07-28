function L = nlap(W)
%generate a normalized graph laplacian from W

for i = 1:size(W,2)
    W(i,i) = 0;
end
d = sqrt(sum(W,2));

for i = 1:size(W,2)
    W(i,:) =W(i,:)/d(i);  
end

for j = 1:size(W,1)
    W(:,j) = W(:,j)/d(j);
end

for i = 1:size(W,2)
    W(i,i) = -1;
end

L = -W;


end