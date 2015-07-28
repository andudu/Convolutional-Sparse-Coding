function [ L ] = ulap( W )
%generate unnormalized graph laplacian
L = -W;
for i= 1:size(L,1)
   W(i,i) = 0;
   L(i,i) = sum(W(i,:),2);   
end


end

