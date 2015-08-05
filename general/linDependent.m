%{
%  linDependent finds linearly independent columns in a given matrix.
%
% X - matrix
% y - vector whos indices are X linearly dependent columns of the first column
%
% example:

A = [1 2 randn;
        2 4 randn;
        3 6 randn]
y = linDependent(A);
A_with_alias_columns_removed = A(:, y)

%
% Dan I. Malta 2013
%}
function y = linDependent(X)
    % first column
    X1 = X(:, 1);
    y = [1];
    
    % iterate thru X columns
    i = 1;
    for j = 1 : size(X, 2)
        % augment next column
        i = i + 1;
        X1 = [X1 X(:, j)];
        
        % test dependence
        if det(X1' * X1) < 1e-6
            X1 = [X1(:, 1 : size(X1, 2) - 1)];
        else
            y = [y, j];
        end
    end
end
