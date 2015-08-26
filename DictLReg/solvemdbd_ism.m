function x = solvemdbd_ism(ah, d,b)

% solvemdbd_ism -- Solve a multiple diagonal block linear system with a 
%                  scaled identity term by iterated application of the 
%                  Sherman-Morrison equation
%
%         The solution is obtained by independently solving a set of linear 
%         systems of the form (see wohlberg-2015-efficient)
%
%                  (rho I + a_0 a_0^H + a_1 a_1^H + ...) x = b
%
%         In this equation inner products and matrix products are taken along
%         the 3rd dimension of the corresponding multi-dimensional arrays; the
%         solutions are independent over the 1st and 2nd (and 4th, if 
%         non-singleton) dimensions.
%   
% Usage:
%       x = solvedbi_sm(ah, rho, b);
%
% Input:
%       ah          Multi-dimensional array containing a^H
%       d           Diagonal d
%       b           Multi-dimensional array containing b
%
% Output:
%       x           Multi-dimensional array containing linear system solution
%
%   
% Author: Brendt Wohlberg <brendt@lanl.gov>  Modified: 2014-12-18
%
% This file is part of the SPORCO library. Details of the copyright
% and user license can be found in the 'Copyright' and 'License' files
% distributed with the library.

a = conj(ah);



K = size(ah,4);
gamma = zeros(size(a), class(a));
delta = zeros([size(a,1) size(a,2) 1 size(a,4)], class(a));


alpha = bsxfun(@rdivide,a(:,:,:,1),d);

%checknan(alpha)

%beta = b/rho;
%d has to be a 3d array
beta = bsxfun(@rdivide, b, d);



clear b;
for k = 1:K,

  gamma(:,:,:,k) = alpha;
  delta(:,:,1,k) =  1 + sum(bsxfun(@times, ah(:,:,:,k), gamma(:,:,:,k)), 3);
  c = sum(ah(:,:,:,k) .* beta, 3);
  d1 = bsxfun(@times, c, gamma(:,:,:,k));
  beta = beta - bsxfun(@rdivide, d1, delta(:,:,1,k));
  if k <= K-1,
    alpha = bsxfun(@rdivide,a(:,:,:,k+1),d);
    for l = 1:k,
      c = sum(ah(:,:,:,l) .* alpha, 3);
      d1 = bsxfun(@times, c, gamma(:,:,:,l));
      alpha = alpha - bsxfun(@rdivide, d1, delta(:,:,1,l));
    end
  end

end

x = beta;

return

end
