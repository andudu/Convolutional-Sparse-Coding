function [ flat_ind ] = rec2flat( rec_ind, imsz, psz )
% translates a flat 1d index to a 2d index by the 
% left->right/ up->down scheme. imsz, psz are scalars for simplicity
% assumes the patch never touches the boundary

if rec_ind(1) < 1
   disp(['flat2rec:','Nonpositive index 1']);
   flat_ind = -1;
   return;
end

if rec_ind(2) < 1
   disp(['flat2rec:','Nonpositive index 2']);
   flat_ind = -1;
   return;
end

ndivc = imsz-psz+1;

if  rec_ind(1) > ndivc
    disp('flat2rec: Index 1 Overflow ');
   flat_ind = -1;
   return;
end


if  rec_ind(2) > ndivc
    disp('flat2rec: Index 2 Overflow ');
   flat_ind = -1;
   return;
end


flat_ind = ndivc*(rec_ind(1)-1)+rec_ind(2);

end

