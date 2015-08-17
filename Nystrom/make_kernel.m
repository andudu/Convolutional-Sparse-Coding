function [kernel] = make_kernel(f)              
 
kernel=zeros(2*f+1,2*f+1);   
for d=1:f    
    value= 1 / (2*d+1)^2 ;    
    for i=-d:d
        for j=-d:d
            kernel(f+1-i,f+1-j)= kernel(f+1-i,f+1-j) + value ;
        end
    end
end
kernel1 = kernel ./ f;

clear kernel

s = 1;
kernel = zeros((2*f+1)^2,1);
for p = 1:2*f+1
    for q = 1:2*f+1
         kernel(s) = kernel1(p,q,:);
         s = s+1;
    end
end