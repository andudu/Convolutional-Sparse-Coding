% -------------------------------------------------------------------------
% yana landa
% yanalanda@gmail.com
%
% computes a feature vector fow each pixel based on patches
% Input:
% input: original image of size m x n x p (p is the number of channels)
% f: the number of pixels in the patch is (2*f+1)^2
% Output:
% v: matrix of neighboring pixels for each pixel in the input image
% -------------------------------------------------------------------------

function v = featureVec(input,f)

[m n p]=size(input); % size of the image
 
v = zeros(m*n,(2*f+1)^2, p); % memory for the output

% replicate the boundaries of the input image via reflction
input2 = padarray(input,f);
 
k = 1;
for i=1:m
    for j=1:n

    i1 = i + f;
    j1 = j + f;

    W = input2(i1-f:i1+f , j1-f:j1+f,:);
    s = 1;
    for p = 1:2*f+1
        for q = 1:2*f+1
            v(k,s,:) = W(p,q,:);
            s = s+1;
        end
    end

    k = k+1;

    end
end

% -------------------------------------------------------------------------
function out = padarray(input, f)
[m, n, p] = size(input);
M = m+2*f;
N = n+2*f;
out = zeros(m+2*f, n+2*f, p);
out(f+1:M-f,f+1:N-f,:) = input;
for j = 1:f
    for i = f+1:M-f
        out(i,j,:) = input(i-f, f-j+1,:);
        out(i,j+f+n,:) = input(i-f, n-j+1,:);
    end    
end

for j = 1:N
    for i = 1:f
        out(i,j,:) = out(2*f-i+1, j,:);
        out(i+f+m,j,:) = out(M-f-i+1, j,:);
    end    
end

        