function [ L ] = winlap( im, ind1, ind2, psz,opt)
%Computes the graph Lacian from image patches starting at ind1 and ending
%at ind2

stpsz = [1,1];
scrop = im(ind1(1):ind2(1)+psz(1)-1, ind1(2):ind2(2)+psz(2)-1);
s1d = imageblocks(scrop,psz,stpsz); 
s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));

if ~isfield(opt, 'SaveMem') %saving memory option
    opt.SaveMem = 0;
end


%%%%%%%%%%%% Computing Graph Laplacian %%%%%%%%%%%%%%%%%%%%%

if ~opt.SaveMem
    
    A = graphgen(s1d, opt.Graph);
    
    if opt.Graph.Laplacian == 'n',
        L = nlap(A);
    end
    if opt.Graph.Laplacian == 'u',
        L = ulap(A);
    end
    
else
    L = savememlapgen(s1d,opt.Graph);
end




end

