%Visualize Graph weights of a single image patch
function [s1,s2] = showweights(rec_ind, opt)

%load info
psz = opt.psz;
stpsz = [1,1];
coefsz = opt.coefsz;
s1 = zeros(coefsz);
imsz = coefsz+psz-[1,1];
flat_ind = rec2flat(rec_ind,imsz,psz,stpsz);

if isfield(opt,'W'), %showimage from given laplacian
    if isfield(opt,'BackgroundImage')
        sb = opt.BackgroundImage;
    else
        sb = ones(imsz);
    end    
    a = opt.W(:,flat_ind);
    s1 = reshape(a,size(s1))';
    s2 = sb(1:coefsz(1),1:coefsz(2)).*s1;    
else
    s = opt.BackgroundImage;
    opt.xind = flat_ind;
    s1d = imageblocks(s,psz,stpsz);
    s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));
    x = s1d(:,flat_ind);
    a = graphgen(s1d,opt,x);     
    s1 = reshape(a,size(s1))';
    s2 = s(1:coefsz(1),1:coefsz(2)).*s1; 
end



end


