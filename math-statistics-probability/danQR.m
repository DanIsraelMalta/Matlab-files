%{
% danQR performs a generalized QR factorization.
%
% A - mxn matrix to be factorized (m >= n >= p)
% B - pxn matrix to be factorized (m >= n >= p) 
%   where:
%   - A = U*L*Q^T
%   - B = S*Q^T
%   - U & Q are orthogonal
%   - L = [0; lower triangular]
%   - S = [0; lower triangular]
% partial - if 'partial' is non zero, then a partial reduction
%           of A is performed - the first 'partial' columns
%           of A are not reduced to triangular form.
%
Dan I. Malta 2009
%}
function [U, Q, L, S] = danQR(A, B, partial)
  % housekeeping
  [m, n] = size(A);
  [p, n1] = size(B);
  if nargin < 3
    partial = 0;
  end
  
  % partial reduction start column
  if partial
    limit = p + 1;
  else
    limit = 1;
  end
  
  % QR factorization
  [Q, S] = qr(B');
  S = S';
  U = eye(m);
  A = A * Q;
  
  % QL factorization
  for i = n : -1 : limit
    % vector-reversal (instead of house-holder)
    temp = A(1 : m - n + i, i);
    temp = temp(end : -1 : 1);
    [v, beta] = gallery('house', temp);
    v = v(end : -1 : 1);
    temp = A(1 : m - n + i, 1 : i);
    A(1 : m - n + i, 1 : i) = temp - beta * v * (v' * temp);
    
    % fill with zeros
    A(1 : m - n + i - 1, i) = zeros(m - n + i - 1, 1);
    
    % upper
    temp = U(:, 1 : m - n + i);
    U(:, 1 : m - n + i) = temp - beta * temp * v * v';
  end
  L = A;
end
