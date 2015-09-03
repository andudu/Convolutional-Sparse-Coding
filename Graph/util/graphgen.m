function [W] = graphgen(s1d, opt, x)
%generate a graph weight matrix from a set of data points s1d

tau = opt.tau;
if nargin < 3
    x = s1d;
end



%%%%%%%%%%%%%%%%%%%%%% Generating Pairwise Distance %%%%%%%%%%%%%%%%%%%%%%%

if strcmp(opt.Metric, 'Euclidean')
    W = sqdist(s1d,x);
    W = exp(-W/(2*tau));
end

if strcmp(opt.Metric, 'Cosine')
    foo = cosdist(s1d,x);
    foo(foo<=0) = 0; %hard thresholding negative values
    W = exp(-abs((1./(foo+.01)-.9901).^(1.3))/tau);
end

if isfield(opt,'Threshold')
    if ~isempty(opt.Threshold)
        W(W<opt.Threshold) = 0;
    end
end


%%%%%%%%%%%%%%%%%%%%%% Generate Graph Weights   %%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    if strcmp(opt.GraphType, 'Full')
        return;
    end
    
    if strcmp(opt.GraphType,'KNearest')
        k = opt.k;
        for iter = 1:size(W,1)
            [v,is] = sort(W(iter,:),'descend');
            W(iter,:) = 0;
            W(iter,is(1:k)) = v(1:k);
        end
        
    end
    
    
    if strcmp(opt.GraphType, 'Window')
        nsz = opt.nsz;
        coefsz = opt.coefsz;
        for iter = 1:size(W,1);
            if nargin == 2
                Ind = patchind(iter,nsz,coefsz);
            else
                Ind = patchind(opt.xind,nsz,coefsz);
            end
            v = W(iter,Ind);
            W(iter,:) = 0;
            W(iter,Ind) = v;
        end
        if size(W,1) == size(W,2)
            W = max(W,W');
        end
        return;
    end
    
    
    if strcmp(opt.GraphType, 'WindowKNearest')
        nsz = opt.nsz;
        coefsz = opt.coefsz;
        k = opt.k;
        for iter = 1:size(W,1);
            Ind = patchind(iter,nsz,coefsz);
            [v,is] = sort(W(iter,Ind),'descend');
            W(iter,:) = 0;
            W(iter,Ind(is(1:k))) = v(1:k);
        end
        if size(W,1) == size(W,2)
            W = max(W,W');
        end
        return;
    end
    
    
    if strcmp(opt.GraphType, 'SpMask') %add a spatial masking (gaussian for now)
        nsz = opt.nsz;
        coefsz = opt.coefsz;
        W_mask = gaussianmask_gen(nsz,coefsz);
        for iter = 1:size(W,1);
            if nargin == 2
                Ind = patchind(iter,nsz,coefsz);
            else
                Ind = patchind(opt.xind,nsz,coefsz);
            end
            v = W(iter,Ind);
            W(iter,:) = 0;
            W(iter,Ind) = v;
        end
        if size(W,1) == size(W,2)
            W = max(W,W');
            W = W_mask.*W;
        else
            W = W_mask(opt.xind,:).*W;
        end
        
        return;
    end
    
    
    if strcmp(opt.GraphType, 'SpInterpolate') %Interpolate a spatial averaging term
        nsz = opt.nsz;
        coefsz = opt.coefsz;
        W_mask = gaussianmask_gen(nsz,coefsz);
        for iter = 1:size(W,1);
            if nargin == 2
                Ind = patchind(iter,nsz,coefsz);
            else
                Ind = patchind(opt.xind,nsz,coefsz);
            end
            v = W(iter,Ind);
            W(iter,:) = 0;
            W(iter,Ind) = v;
        end
        if size(W,1) == size(W,2)
            W = max(W,W');
            W = W_mask.*W;
        else
            W = W_mask(opt.xind,:).*W;
        end
        
        return;
    end
    











