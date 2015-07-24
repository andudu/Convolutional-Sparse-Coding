function [ flat_ind ] = rec2flat(rec_ind, imsz, psz, stpsz)
% translates a flat 1d index to a 2d index by the 
% left->right/ up->down scheme. 
% assumes the patch never touches the boundary
% rec_ind: rectangular coordinates: (i,j)
% imsz: image size: (m x n)
% psz: patch size: (a x b)
% stpsz: step size: (st1 x st2)
% all indices are row first col second

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

n = imsz(2);
m = imsz(1);
a = psz(1);
b = psz(2);

if  rec_ind(1) > m-psz+1
   disp('flat2rec: Index 1 Overflow ');
   flat_ind = -1;
   return;
end


if  rec_ind(2) > n-psz+1
    disp('flat2rec: Index 2 Overflow ');
   flat_ind = -1;
   return;
end


npc = floor((n-b)/stpsz(2))+1;
nr = floor((rec_ind(1)-1)/stpsz(1));
nc = floor((rec_ind(2)-1)/stpsz(2))+1;

flat_ind = npc*nr+nc;

end

