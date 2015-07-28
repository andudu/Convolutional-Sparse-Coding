function [L,scrop] = graphgen(sh, wsz, psz,neig,tau)

%generate a graph laplacian for each patch with sliding windows

%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%% Divide the Windows %%%%%%%%%%%%%%%%
L = {};  
ind1 = [1,1];
ind2 = ind1+wsz-[1,1];
k = 1;

while(ind2(1)<= size(sh,1))
    while(ind2(2)<= size(sh,2))
        L{k}.ind1 = ind1;
        L{k}.ind2 = ind2;
        a = winlap(sh,ind1,ind2,psz,tau);
        [phi,e] = eigs(a,neig,'sr');
        e = diag(e);
        L{k}.phi = phi;
        L{k}.E = e;
        ind1(2) = ind1(2)+wsz(2);
        ind2 = ind1+wsz-[1,1];  
        k = k+1;
    end
    ind1(2)= 1;
    ind1(1) = ind1(1) + wsz(1);
    ind2 = ind1+wsz-[1,1]; 
end

% %%%%%%%%%%%%%%% Post Processing %%%%%%%%%%%%%%%%%%%%%%
x_max = L{end}.ind2(1);
y_max = L{end}.ind2(2);
scrop = sh(1:x_max,1:y_max);
% sh = scrop;
% sl = sl(1:x_max,1:y_max);
% 
% tag = 'noise';
% save(['data_graph_',tag,'.mat'],'L','sh','sl');
% 
% 
% x = [];
% for i = 1:16
%     temp = L{i}.phi(:,30);
%     l1 = L{i}.ind1(1):L{i}.ind2(1);
%     r1 = L{i}.ind1(2):L{i}.ind2(2);
%     x(:,:,i) = reshape(temp,length(l1),length(r1));
% end
% square_plot(x,{});