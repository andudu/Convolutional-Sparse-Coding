%generate eigenvectors for the flicker dataset via Nystrom;

load('Flicker1_512_split.mat');
imnum = size(S,3);
lambda = 4;


for i = 1:imnum
    i
    temp = double(S(:,:,i))/255;
    [~,sh] = lowpass(temp,lambda,15);
    clear temp;
    imsz = size(sh);
    psz = [12,12];
    stpsz = [1,1];
    sh = padarray(sh,psz-[1,1],'symmetric','post');
    
    %generate graph
    opt.tau = 1.2;
    opt.Laplacian = 'n';
    opt.numsample = 800;
    opt.Metric = 'Euclidean';
    opt.neig = 30;
    data = imageblocks(sh,psz,stpsz);
    data = reshape(data,size(data,1)*size(data,2),size(data,3));
    data = data';
    [phi,E] = nystrom(data,opt);
    save(['Eigenvectors/eig',num2str(i),'.mat'],'phi','E');
    
end

imsz = [size(S,1),size(S,2)];
save('Eigenvectors/param.mat','imnum','lambda','imsz','psz');
    

    