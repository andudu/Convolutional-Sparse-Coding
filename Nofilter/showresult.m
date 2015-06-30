function showresult(X,D,opt)
%quick way of displaying useful info


scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

%plotting full reconstruction                           
figure;

DX = scnv(D, X);
if isfield(opt,'lowfreq')
    sl = opt.lowfreq;
    DX = DX+sl;
end

imagesc(DX);
colormap(gray);
axis off;
title ('full reconstruction');

%plotting nonzero elements
plot_opt.sparse = 1;
square_plot(X,plot_opt);

%plotting the reconstruction componentwise
temp = zeros(size(X));
for i = 1:size(X,3)
    temp(:,:,i) = ifft2(bsxfun(@times,fft2(D(:,:,i),size(X,1),size(X,2)),...
        fft2(X(:,:,i))));
end
plot_opt.sparse = 0;
plot_opt.grey = 1;
square_plot(temp,plot_opt);


%plotting the distribution of coefficient.
coeff = reshape(sum(sum(abs(X),1),2),1,size(X,3));
figure;
plot(1:1:size(X,3),log(coeff));
title('distribution of log coeff sum');




