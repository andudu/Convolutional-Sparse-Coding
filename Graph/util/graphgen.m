function [W] = graphgen(s1d, opt)
%generate a graph weight matrix from a set of data points s1d

tau = opt.tau;

%%%%%%%%%%%%%%%%%%%%%% Generating Pairwise Distance %%%%%%%%%%%%%%%%%%%%%%%
if strcmp(opt.Metric, 'Euclidean')
    W = sqdist(s1d,s1d);
    W = exp(-W/(2*tau)); 
end

if strcmp(opt.Metric, 'Cosine')
    foo = cosdist(s1d,s1d);
    foo(foo<-0) = 0; %hard thresholding negative values
    W = exp(-(1./(foo+.01)-.9901)/tau);
end


%%%%%%%%%%%%%%%%%%%%%% Generate Graph Weights   %%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(opt.GraphType, 'Full')
    return;
end

if strcmp(opt.GraphType, 'Window')
    nsz = opt.nsz;
    coefsz = opt.coefsz;
    for iter = 1:size(W,1);
        Ind = patchind(iter,nsz,coefsz);
        v = W(Ind,iter);
        W(:,iter) = 0;
        W(Ind,iter) = v;
    end
    W = max(W,W');  
    return;
end


if strcmp(opt.GraphType, 'WindowKNearest')
    nsz = opt.nsz;
    coefsz = opt.coefsz; 
    k = opt.k;
    for iter = 1:size(W,1);
        Ind = patchind(iter,nsz,coefsz);
        [v,is] = sort(W(Ind,iter),'descend');
        W(:,iter) = 0;
        W(Ind(is(1:k)),iter) = v(1:k);
    end
    W = max(W,W');      
end


if strcmp(opt.GraphType, 'SoftWindow')
    if ~isfield(opt,'SpatialWeight')
        Ws = opt.SpatialWeight ;
    else
        [a,b] = meshgrid(1:1:coefsz(1),1:1:coefsz(2));
        c = reshape(cat(3,a,b),coefsz(1)*coefsz(2),2)';
        Ws = sqdist(c,c);
        Ws = exp(-Ws/10);
    end
        W = Ws.*W;
    return; 
end


