%{
% invDet returns the power (of a given order) of the inverse determinant of of a symmetric (real or complex) positive definite matrix.
%
% xi_mat - symmetric (real or complex) positive definite matrix
% xi_power - required power (must be positive real)
% xo_inv - |{xi_mat}|.^{xi_power}. if matrix is not symmetric & PD - return Inf.
%
% Dan I. Malta (2016)
%}
function xo_inv = invDet(xi_mat, xi_power)
    % locals
    n = size(xi_mat, 1);
   
    % validity check
    if ~ismatrix(xi_mat)
        error('invDet: xi_mat must be a matrix.');
    end
    if ndims(xi_mat) > 2 %#ok<ISMAT>
        error('invDet: xi_mat must be a 2D matrix.');
    end
    if size(xi_mat, 2) ~= n
        error('invDet: xi_mat must be square.');
    end
    if nargin == 1
        xi_power = 1;
    end
    if ~isnumeric(xi_power) || ~isreal(xi_power) || (numel(xi_power ) ~=  1) || (xi_power <= 0)
        error('invDet: xi_power must be positive real number.');
    end
    % does matrix includes elements which are zero?
    if nnz(xi_mat - xi_mat') ~= 0
        xo_inv = Inf;
    else
        % cholsky decomposition
        [R, q] = chol(xi_mat);
       
        % is matrix positive definite?
        if q == 0
            xo_inv = prod(diag(R)).^(-xi_power);
        else % Cholsky say's matrix is not positivie definite
            % check using eigenvalues
            eigs = eig(xi_mat);
            if any(eigs <= 0)
                xo_inv = Inf;
            else
                xo_inv = prod(eigs).^(-xi_power);
            end
        end
    end
end
