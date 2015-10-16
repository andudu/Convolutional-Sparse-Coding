% generate and format for time test


% load the eigenvectors. 



load(strcat(sporco_path,'/Graph/CacheData/Lena_Knearest_vark/Mat/EuclideanM2.mat'));
temp = L; 
clear L; 
A = temp{1,1}.M; 
[V,E] = eigs(A,400,'sr'); 
E = diag(E); 
mkdir(strcat(sporco_path,'/Graph/CacheData/Lena_timetest')); 
for i = 1:4
    L = {}; 
    L{1}.ind1 = temp{1,1}.ind1;
    L{1}.ind2 = temp{1,1}.ind2;
    numeig = 100*i; 
    L{1}.phi = V(:,1:1:numeig); 
    L{1}.E = E(1:1:numeig); 
    save(strcat(sporco_path,'/Graph/CacheData/Lena_timetest/Eig',num2str(i),'.mat'),'L','numeig'); 
end


% L2{1}.ind1 = L{1,1}.ind1; 
% L2{1}.ind2 = L{1,1}.ind2; 
% L1{1}.phi = V; 
% L1{1}.E = E; 
% L2{1}.M = A; 
% clear A V E L; 



% load the full matrix. 


