function square_plot(D,opt)

% D is m*n*p size dictionary
% opt is options
% opt.grey for greyscale
% opt.sparse plots nonzero elements

p = size(D,3);

k = ceil(sqrt(p));

figure;

if nargin == 0 ,
    opt = {};
end

if~isfield(opt,'sparse'),
    opt.sparse = 0;
    if ~isfield(opt,'grey'),
        opt.grey = 0;
    end
end



for i = 0:k-1
    for j = 1:k
        ind = i*k+j;
        if ind <= p,
            subplot(k,k,ind);
            if ~opt.sparse
                if ~opt.grey
                    imagesc(D(:,:,ind));
                    axis off;
                else
                    imagesc(D(:,:,ind));
                    colormap(gray);
                    axis off;
                end
            else
                spy(D(:,:,ind));
                axis off;
            end
        end
    end
end

