function [ rec_ind ] = flat2rec( flat_ind, imsz, psz, stpsz )
% translates a flat 1d index to a 2d index by the 
% left->right/ up->down scheme. 
% assumes the patch never touches the boundary
% rec_ind: rectangular coordinates: (i,j)
% imsz: image size: (m x n)
% psz: patch size: (a x b)
% stpsz: step size: (st1 x st2)
% all indices are row first col second

if flat_ind < 1
   disp(['flat2rec:','Nonpositive index']);
   rec_ind = -1;
   return;
end

n = imsz(2);
m = imsz(1);
a = psz(1);
b = psz(2);

npc = floor((n-b)/stpsz(2))+1;
npr = floor((m-a)/stpsz(1))+1;

if flat_ind > npc*npr
   disp('flat2rec: Index Overflow ');
   rec_ind = -1;
   return;
end



nr = floor((flat_ind-1)/npc); %number of rows above current position
nc = flat_ind-nr*npc;   %number of patches in the last row til this patch

rec_ind(1) = nr*stpsz(1)+1;
rec_ind(2) = (nc-1)*stpsz(2)+1;


