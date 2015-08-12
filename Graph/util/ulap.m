function [ L ] = ulap( W )
%generate unnormalized graph laplacian
L = -W;
avgd = 0;
for i= 1:size(L,1)
   W(i,i) = 0;
   foo = sum(W(i,:),2);   
   L(i,i) = foo;
   avgd = avgd + foo;
end
avgd = avgd/size(L,1);
L = L/avgd;

end

