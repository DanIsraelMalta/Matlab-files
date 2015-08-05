%{
% danChol performs the Cholesky factorization for positive semidefinite matrices
%
% danChol can handle semi positive matrix, while matlab builtin chol requires positive definite matrix.
%
% A - symmetric positive semi-definite matrix
% L - lower triangular cholesky factor matrix (A = L * L')
% def -  -1 = A was negative definite
%               0 = A was positive semi definite
%               1 = A was positive definite
%
% Dan I. Malta 2014
%}
function [L, def] = danChol(A)
    % housekeeping
    L  = zeros(size(A));
    def = 1;
    
    % iterate over rows
    for i = 1 : size(A, 1)
        for j = 1 : i
            s = A(i,j);
            for k = 1 : j - 1
                s = s - L(i, k) * L(j, k);
            end
            if j < i % A is positive semi definite
                if L(j, j) > eps
                    L(i,j) = s / L(j, j);
                else
                    L(i, j) = 0;
                end
            else % A is negative definitie
                if s < -eps
                    s = 0;
                    def = -1;
                elseif s < eps
                    s = 0;
                    def = min(0, def);
                end
                L(j,j) = sqrt(s);
            end
        end
    end
end
