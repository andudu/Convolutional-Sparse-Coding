function [ Y ,rate] = compress_nd( X )
%Takes a 2d or 3d or 4d array and compress into a cell of sparse matrices
%converts singles to doubles
%rate is the percent of nonzero entries.

xsz = size(X);
rate =0;

ndim = length(xsz);
if ndim == 2,
    Y = {};
    Y{1} = sparse(double(X));
    rate = rate + nnz(Y{1});
end
if ndim == 3,
    Y = cell(xsz(3),1);
    for i = 1:xsz(3)
        Y{i} = sparse(double(X(:,:,i)));
            rate = rate + nnz(Y{i});
    end
end
if ndim == 4,
    Y = cell(xsz(3),xsz(4));
    for i = 1:xsz(3)
        for j = 1:xsz(4)
            Y{i,j} = sparse(double(X(:,:,i,j)));
                rate = rate + nnz(Y{i,j});
        end
    end
end
rate = rate/numel(X);




        
        



end

