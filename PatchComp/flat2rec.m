function [ rec_ind ] = flat2rec( flat_ind, imsz, psz )
% translates a flat 1d index to a 2d index by the 
% left->right/ up->down scheme. imsz, psz are scalars for simplicity
% assumes the patch never touches the boundary

if flat_ind < 1
   disp(['flat2rec:','Nonpositive index']);
   rec_ind = -1;
   return;
end

ndivc = imsz-psz+1;

if flat_ind > ndivc*ndivc
    disp('flat2rec: Index Overflow ');
   rec_ind = -1;
   return;
end


rec_ind(1) = floor((flat_ind-1)/ndivc)+1;
rec_ind(2) = flat_ind - (rec_ind(1)-1)*ndivc;


end

