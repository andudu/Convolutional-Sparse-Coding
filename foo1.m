function [D] = foo1(X, Y)
%usage: 
% X = (x1,x2,\dots, xm) m col vectors
% Y = (y1,y2,\dots, yn) n col vectors
% D_ij = |xi-yj|^2 

m = size(X,2);
n = size(Y,2);  
Yt = Y';
XX = sum(X.*X,1);        
YY = sum(Y.*Y,1);
YYt = YY';
D = XX(ones(1,n),:) + YYt(:,ones(1,m)) - 2*Yt*X;
end

