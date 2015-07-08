%timetest for cropping

z = rand(1000,1000);
y = [];

tic;
for i = 1:3000
    y = z(1:500,1:500);
    y = padarray(y,[500,500],'post');
end
toc;


tic;
z_times = ones(500,500);
z_times = padarray(z_times,[500,500],'post');
for i = 1:1000
    y=z;
    y = y.*z_times;
end
toc;
