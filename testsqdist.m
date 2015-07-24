%test the graphs

%%%%%%%%%%%%%%  Loading Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s = single(rgbtogrey(stdimage('lena')))/255;
% s = imresize(s,.5);
% 
% s = s+.1*randn(size(s));

load('noise_data.mat');

[sl,sh] = lowpass(s,5,15);

%%%%%%%%%%%%  Setting Patch Info %%%%%%%%%%%%%%%%%%%%%%%%%%
psz = [12,12];
stpsz = [1,1];
scrop = sh(200:256,110:165);
imsz = size(scrop);
s1d = imageblocks(scrop,psz,stpsz);
coefsz = imsz-psz+[1,1];
s1d = reshape(s1d,size(s1d,1)*size(s1d,2),size(s1d,3));


%%%%%%%%%%%% Computing Graph Laplacian %%%%%%%%%%%%%%%%%%%%%
A = sqdist(s1d,s1d);
numnb = 4;
%Add a mask, each patch only compares with 8 neighbors
% zcrop = zeros(size(A));
% for i = 1:size(A,1)
%     i0 = flat2rec(i,imsz,psz,stpsz);
%     nb_ind = pnbr(i0,numnb,coefsz);
%     for j = 1:size(nb_ind,2)
%         itemp = rec2flat(nb_ind(:,j),imsz,psz,stpsz);
%         zcrop(i,itemp) = 1;
%         zcrop(itemp,i) = 1;
%     end
% end    
% 
% A(A>0.15) = 0;
A = exp(-A/2);
% A = A.*zcrop;
for i = 1:size(A,2)
     A(i,i) = 0;
%     A(i,:) = A(i,:) / sum(A(i,:));
%     A(i,i) = -1;
    A(i,i) = -sum(A(i,:));
end
L = -A;



%%%%%%%%%%%%%%%% Computing and Plotting Eigenvalues %%%%%%%%
eignum = 10;
[v,e] = eigs(L,eignum,'sr');


for i = 1:eignum
    temp = zeros(coefsz(1),coefsz(2));
    for j = 1:size(s1d,2) 
        ind =flat2rec(j,imsz,psz,stpsz); 
        temp(ind(1),ind(2)) = v(j,i);
    end
    figure;
    imagesc(temp);
    colorbar;
    title([num2str(i),'th eig']);
end



