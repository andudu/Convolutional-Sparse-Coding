%compute the cosine distance betwen patches

% Construct initial dictionary
% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
D = D(:,:,6); %get the 6th element
D = double(D);
load ('noise_data.mat'); %loading the noise
s = double(s);

a = [];
for i = 1:size(s,1)-12
    for j = 1:size(s,2)-12
        temp = s(i:i+11,j:j+11);
        a(i,j) = (temp(:))'*D(:);
    end
end
a = abs(a);
a(a<0.2)= 0 ;

figure;
imagesc(a);
colorbar;

