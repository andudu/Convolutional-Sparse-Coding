crpsz = 256;

image_names = { 'airplane.png' ,   'lena.grey.png', 'boats.png' ,'house.png' ,...
    'monarch.png','barbara.grey.png',  'bridge.grey.png',   'kiel.grey.png',...
    'man.grey.png', 'peppers.png', 'goldhill.png',...
    'mandrill.png' };
%image_num = size(image_names,2);
%image_index =  [2,6,7,8,9,10,11];

image_index = 1; %barbara
image_num = length(image_index);

S0 = zeros(crpsz, crpsz, image_num , 'single');

ind = 1;
for i = image_index
    tmp = stdimage(image_names{i});
    if ndims(tmp) == 3
        tmp = rgb2gray(tmp);
    end
    tmp = imresize(tmp,[crpsz,crpsz]);
    S0(:,:,ind) = single(tmp)/255;
    ind = ind+1;
end
clear tmp tmp1; 

npd = 16;
Sh = zeros([size(S0),8]);
Sl = zeros([size(S0),8]);
diff = zeros([size(S0),7]);
for i = 1:8
    lambda = 1*i;
    [Sl(:,:,i), Sh(:,:,i)] = lowpass(S0, lambda, npd);
    if i < 8,
        diff(:,:,i) = Sh(:,:,i+1)-Sh(:,:,i);
    end
end

o1.grey = 1;
o1.unifscale = 0;
square_plot(Sh,o1);
square_plot(Sl,o1);
square_plot(diff,o1);