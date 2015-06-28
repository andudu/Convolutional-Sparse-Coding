function X = decompress_nd(Y)
%decompress cells of sparse matrices into corresponding 4d or 3d arrays

ysz = size(Y);
ndim = length(ysz);
xsz = size(Y{1});

if ndim == 2,
    xsz = [xsz,ysz];
    X = zeros(xsz);
    for i =1:ysz(1)
        for j = 1:ysz(2)
            X(:,:,i,j) = single(full(Y{i,j}));
        end
    end
end

if ndim == 1,
    if numel(ysz) == 1 
        X = single(full(Y{i,j}));
    else
        xsz = [xsz, numel(ysz)];
        X = zeros(xsz);
        for i = 1:numel(ysz)
            X(:,:,i) = single(full(Y{i}));
        end
    end
end

            
            