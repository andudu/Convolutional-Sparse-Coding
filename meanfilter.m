function y = meanfilter(s,stpsz,psz)
%script for computing mean of each patch and averaging 
% s is the 2d signal with dimensions nx x ny
% y is the patch average of the block mean
% stpsz is the stepsize taken right and down for each patch
% psz is the patch size
%(not taking care of boundary cases yet)

nx = size(s,1);
ny = size(s,2);
s_div = zeros(size(s));
y = zeros(size(s));
d = psz*psz;


for i = 1:stpsz:nx-stpsz+1
        i_end = min([i+psz-1,nx]);
    for j = 1:stpsz:ny-stpsz+1
        j_end = min([j+psz-1,ny]);
        patch_avg = sum(vec(s(i:i_end,j:j_end)))/d;
        s_div(i:i_end,j:j_end) = s_div(i:i_end,j:j_end) + ones(i_end-i+1,j_end-j+1);
        y(i:i_end,j:j_end) = y(i:i_end,j:j_end) + patch_avg;
    end
end
y = y./s_div;



