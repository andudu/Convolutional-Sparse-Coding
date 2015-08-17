% -------------------------------------------------------------------------
% yana landa
% yanalanda@gmail.com
%
% Input: 
% data: feature vector (matrix), the third dimension is the number of
% channels in the image (should be modified if use different type of data)
% kernel: gaussian to convolve with the data
% num_samples: number of eigenvalues/vectors, also the number of random
% samples
% sigma: value used in computing similarity
% Output: 
% V: eigenvectors
% L: eigenvalues of Laplacian
%
% -------------------------------------------------------------------------

function [V L] = nystrom_color(data, kernel, num_samples, sigma)

% randomly select samples
num_rows = size(data, 1);
permed_index = randperm(num_rows);
sample_data = data(permed_index(1:num_samples), :, :);
other_data = data(permed_index(num_samples+1:num_rows), :, :);
clear data;

% calculate the euclidean distance between samples themselves
disp('Calculating A');
A = zeros(num_samples, num_samples);
sigma = sigma*sigma;
for i = 1:num_samples
    for j = 1:num_samples
        d1 = (sample_data(i,:,:)-sample_data(j,:,:)).^2;
        d2 = sum(d1,3);
        d = sum(kernel'.*d2);
        A(i,j) = exp(-d/sigma);
    end
end
A = single(A);

% calculate the euclidean distance between samples and other points
disp('Calculating B');
other_points = num_rows - num_samples;
B = zeros(num_samples, other_points);
for i = 1:num_samples
    for j = 1:other_points
        d1 = (sample_data(i,:,:) - other_data(j,:,:)).^2;
        d2 = sum(d1,3);
        d = sum(kernel'.*d2);

        B(i,j) = exp(-d/sigma);
    end
end
B = single(B);
clear sample_data other_data;


% Normalize A and B using row sums of W, where W = [A B; B' B'*A^-1*B].
% Let d1 = [A B]*1, d2 = [B' B'*A^-1*B]*1, dhat = sqrt(1./[d1; d2]).

disp('Normalizing A and B for Laplacian...');
B_T = B';
d1 = sum(A, 2) + sum(B, 2);
d2 = sum(B_T, 2) + B_T*(pinv(A)*sum(B, 2));
dhat = sqrt(1./[d1; d2]);
A = A .* (dhat(1:num_samples)*dhat(1:num_samples)');
B1 = dhat(1:num_samples)*dhat(num_samples+(1:other_points))';
B = B .* B1;
clear d1 d2 B1 dhat;

% Do orthogalization and eigendecomposition
disp('Orthogalizing and eigendecomposition...');
Asi = sqrtm(pinv(A));
B_T = B';
BBT = B*B_T;
W = single(zeros(size(A, 1)+size(B_T, 1), size(A, 2)));
W(1:size(A, 1), :) = A;
W(size(A, 1)+1:size(W, 1), :) = B_T;
clear B B_T;
% Calculate R = A + A^-1/2*B*B'*A^-1/2
R = A + Asi*BBT*Asi;
R = (R + R')/2; % Make sure R is symmetric, sometimes R can be non-symmetric because of numerical inaccuracy
[U L] = eig(R);
[val ind] = sort(diag(L), 'descend');
U = U(:, ind); % in decreasing order
L = L(ind, ind); % in decreasing order
clear A R BBT;
W = W*Asi;
V = W*U(:, 1:num_samples)*pinv(sqrt(L(1:num_samples, 1:num_samples)));

V(permed_index,:) = V;
V = real(V);
L = 1-diag(L);
