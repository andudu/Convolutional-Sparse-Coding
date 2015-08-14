%Visualize Graph weights of a single image patch 
function [s1,s2] = showweights(rec_ind, W, opt)

%load info
imsz = opt.imsz;
psz = opt.psz;
stpsz = [1,1];
if isfield(opt,'BackgroundImage')
    sb = opt.BackgroundImage;
else 
    sb = ones(imsz);
end
coefsz = imsz - psz + [1,1];
s1 = zeros(coefsz);

flat_ind = rec2flat(rec_ind,imsz,psz,stpsz);
a = W(:,flat_ind);
s1 = reshape(a,size(s1))';
s2 = sb(1:coefsz(1),1:coefsz(2)).*s1;

