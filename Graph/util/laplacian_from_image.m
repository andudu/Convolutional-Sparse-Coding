function [L,scrop] = laplacian_from_image(im,opt)
%function for generating a series of windowed graph laplacian
%from a givin image. L is formatted for the cbpdn_L code.  
%Format for L: L{i}.ind1, L{i}.ind2 L{i}.M (or L{i}.phi, L{i}.E)

%Format for opt:
%   opt.Lformat: 'Full', 'Sparse', 'Eig'
%   opt.neig: number of eigenvectors to save if 'Eig' option
%   opt.Laplacian: 'n', 'u'
%   opt.Graph:
%           opt.Graph.tau
%           opt.Graph.Metric: 'Euclidean', 'Cosine'
%           opt.Graph.GraphType: 'Full', 'Window', 'WindowKNearest',
%                                'SoftWindow'
%           opt.Graph.nsz
%           opt.Graph.coefsz (set by default here, from now on, equal to imsz)
%           opt.Graph.k: number of k-nearest neighbors


%pad the image first
im = padarray(im,opt.psz-[1,1],'symmetric','post');
if ~isfield(opt,'Nystrom')
    opt.Nystrom = 0;
end

imsz = size(im);
wsz = opt.wsz;
psz = opt.psz;
opt.Graph.coefsz = wsz;

n1 =  floor((imsz(1)+1-psz(1))/wsz(1)); %number of windows in a col
n2 =  floor((imsz(2)+1-psz(2))/wsz(2)); %number of windows in a row
coefsz = [n1*wsz(1), n2*wsz(2)];

iter = 1;
L = {};


if strcmp(opt.Graph.GraphType, 'SoftWindow')
    [a,b] = meshgrid(1:1:wsz(1),1:1:wsz(2));
    c = reshape(cat(3,a,b),wsz(1)*wsz(2),2)';
    Ws = sqdist(c,c);
    opt.Graph.SpatialWeight = exp(-Ws/10);
end


for i = 1:n1
    for j = 1:n2
        ind1 = [(i-1)*wsz(1)+1,(j-1)*wsz(2)+1];
        ind2 = ind1 + wsz - [1,1];
        L{iter}.ind1 = ind1;
        L{iter}.ind2 = ind2;
        
        if opt.Nystrom  %skip graph generation and use nystrom to do that
            stpsz = [1,1];
            scrop = im(ind1(1):ind2(1)+psz(1)-1, ind1(2):ind2(2)+psz(2)-1);
            s1d = imageblocks(scrop,psz,stpsz);
            s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));
            o.tau = o.Graph.tau; o.Metric = o.Graph.Metric;
            o.Laplacian = opt.Laplacian; o.neig = opt.neig; 
            o. numsample = 500; %tentative number of samples
            [V,E] = nystrom(s1d',opt);    
                L{iter}.E = diag(E);
                L{iter}.phi = V;            
        else
            Ltemp = winlap( im, ind1, ind2, psz,opt);
            if strcmp(opt.Lformat, 'Full')
                L{iter}.M = Ltemp;
            end
            if strcmp(opt.Lformat, 'Sparse')
                L{iter}.M = sparse(Ltemp);
            end
            if strcmp(opt.Lformat, 'Eig')
                [V,E] = eigs(Ltemp,opt.neig,'sr');
                L{iter}.E = diag(E);
                L{iter}.phi = V;
            end
        end
        iter = iter+1;
    end
end


%cropped image that fits this windowing
scrop = im(1:coefsz(1), 1:coefsz(2));

end
