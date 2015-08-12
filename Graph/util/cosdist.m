function [D] = cosdist(X, Y)
%usage: 
% X = (x1,x2,\dots, xm) m col vectors
% Y = (y1,y2,\dots, yn) n col vectors
% D_ij = |xi-yj|^2 
 
Yt = Y';  
XX = sqrt(sum(X.*X,1));        
YY = sqrt(sum(Yt.*Yt,2));   
XX(XX == 0) = .01;
YY(YY ==0 ) = .01; %division by 0
D = bsxfun(@rdivide,Yt,YY) * bsxfun(@rdivide,X,XX);

end

