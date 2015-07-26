function [ L ] = foowinlap( im, ind1, ind2, psz )
%computes the graph laplacian for a given window
%does not do boundary checks

stpsz = [1,1];
scrop = im(ind1(1):ind2(1)+psz(1)-1, ind1(2):ind2(2)+psz(2)-1);
s1d = imageblocks2(scrop,psz,stpsz); %down then right for compatibility with reshape

%add a distance term?

%should I use normalized or unnormalized? 

s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));


%%%%%%%%%%%% Computing Graph Laplacian %%%%%%%%%%%%%%%%%%%%%
A = sqdist(s1d,s1d);
A = exp(-A/2);
for i = 1:size(A,2)
    A(i,i) = 0;
end
d = sqrt(sum(A,2));

for i = 1:size(A,2)
    A(i,:) =A(i,:)/d(i);  
end

for j = 1:size(A,1)
    A(:,j) = A(:,j)/d(j);
end

for i = 1:size(A,2)
    A(i,i) = -1;
end

L = -A;


end

