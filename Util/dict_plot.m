function h = dict_plot(D,opt)
% D is m*n*p size dictionary
% opt is options
% opt.grey for greyscale
% opt.unifscale for uniform scale
% opt.sz for cropping

p = size(D,3);
k = ceil(sqrt(p));

figure;
if nargin == 1 ,
    opt.foo = 0;
end


if ~isfield(opt,'grey'),
    opt.grey = 1;
end
if ~isfield(opt,'unifscale'),
    opt.unifscale = 1;
end
if ~isfield(opt,'fsz'),
    opt.fsz = repmat([size(D,1);size(D,2)],[1,size(D,3)]);
end
if ~isfield(opt,'axis'),
    opt.axis = 1;
end


if opt.unifscale,
    m1 = max(vec(D));
    m2 = min(vec(D));
    D = (D-m2)*64/(m1-m2);
end

for i = 0:k-1
    for j = 1:k
        ind = i*k+j;
        if ind <= p,
            subplot(k,k,ind);
            if ~opt.axis
                axis off;
            end
            temp = D(:,:,ind);
            temp = temp(1:opt.fsz(1,ind),1:opt.fsz(2,ind));
            if opt.unifscale,
                image(temp);
            else
                imagesc(temp);
            end                
            if opt.grey
                colormap(gray);
            end
        end
    end
end

h = gcf;