%{
% orthogonalPolyFit fit's an orthogonal polynomial to a time series
%
% I find orthogonal polynomials to be better then the linear least square
% method employed by matlab's polyfit, since they do no involve any matrix
% linear system solution, thuse they:
% 1) don't suffer from numeric instability
% 2) can be applied recursively and efficiently for large data sets (perfect for smoothing purposes).
% 3) can use probabilistic criteria's to ensure their convergence.
%
% x, y - time sereis {x, y}
% n      - polynom degree
% ys   - orthogonal polynom values at given x
% p     - orthogonal polynom coefficients
%
% example:

x = linspace(0, 10, 300);
y = sin(x.^3 / 100) .^ 2 + 0.05 * randn(size(x));
n = optimalPoly(x,y);
ys = orthogonalPolyFit(x, y, n);
plot(x, y, '.', x, ys, 'k');
grid on;
title(['optimal degree of polynom is ', num2str(n)]);

%
% Dan I. Malta 2009
%}
function [ys, p] = orthogonalPolyFit(x, y, n)
    % housekeeping
    p = mean(y(:));
    ys = ones(size(y)) * p;

    % columnize
    x = x(:);
    siz0 = size(y);
    y = y(:);

    % orthogonal polynomials coefficents
    p = zeros(3, n + 1);
    p(2, 2) = mean(x);
    N = length(x);
    PL = ones(N, n + 1);
    PL(:, 2) = x - p(2, 2);
    for i = 3 : n + 1
        p(2, i) = sum(x .* PL(:, i - 1).^2) / sum(PL(:, i - 1).^2);
        p(3, i) = sum(x .* PL(:, i - 2) .* PL(:, i - 1)) / sum(PL(:, i - 2).^2);
        PL(:, i) = (x - p(2, i)) .* PL(:, i - 1) - p(3, i) * PL(:, i - 2);
    end
    p(1, :) = sum(repmat(y, [1, n + 1]) .* PL) ./ sum(PL.^2);
     
    % smooth output
    ys = sum(PL .* repmat(p(1, :), [N, 1]), 2);
    ys = reshape(ys, siz0);
     
    % Coefficients of the polynomial in its final form
    yi = zeros(n + 1);
    yi(1, n + 1) = 1;
    yi(2, n : n + 1) = [1 -p(2,2)];
    for i = 3 : n + 1
        yi(i, :) = [yi(i - 1, 2 : end) 0] - p(2, i) .* yi(i - 1, :) - p(3, i) * yi(i - 2, :);
    end
    p = sum(repmat(p(1, :).', [1, n + 1])  .*  yi, 1);
end
