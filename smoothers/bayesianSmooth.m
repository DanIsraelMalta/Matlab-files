%{
% bayesianSmooth minimizes a quadratic cost function by utilizing both the
% maximum likelihhos and the apriori pdf (a complete bayesian estimator) - thus returning a perfect
% aposteriori pdf! this should, in theory and for non linear observations,
% perform much better then any optimal-linear or  maxmim-likelihhood based smoother.
%
% in theory, this smoother should give good results in environments with signal-to-noise ratio of 0.3 (equivaent to -10db).
%
% bayesianSmooth minimizes the same cost function like directSmooth.
%
% x - signal to be smoothed
% order - signal derivative order to be smoothed (2 or 3)
% w - weight given to the derivative smoothing component in the costfunction
% y - smoothed output
%
% example:

t = 0 : 0.01 : 10;
x = cos(t);
y = x + randn(size(t));
z = bayesianSmooth(y, 3, 1e9);
snrDB = 20*log10((norm(x) / norm(y))^2);

figure;
plot(t, y, 'r', t, x, 'k', t, z, 'b');
legend('observation', 'original signal', 'smoothed observation')
title(['signal-to-noise ratio is ', num2str(snrDB), '[db]']);
grid on;

%
% Dan I. Malta 2015
%}
function y = bayesianSmooth(x, order, w)
    % columnize
    x=x(:);
    
    % housekeeping
    n = length(x);
    A0 = eye(n); % part corresponding to distance function

    % writing cost function in quadratic matrix form: min(yAy + yb + c)
    switch order
        case 2 % second derivative component
            A1 = diag([1 5 6 * ones(1, n - 4) 5 1]) - 4 * diag([1 2 * ones(1, n - 3) 1], 1) + 2 * diag(ones(1, n - 2), 2);
        otherwise % third derivative component (any higher order demand, shall be treated as third order)
            A1 = diag([1 10 19 20 * ones(1, n - 6) 19 10 1]) - 6 * diag([1 4 5 * ones(1, n - 5) 4 1], 1) +...
                      6 * diag([1 2 * ones(1, n - 4) 1], 2) - 2 * diag(ones(1, n - 3), 3);
    end
    
    % combine distance and derivative components into a symmetric matrix
    A = A0 + w * A1;
    A = 0.5 * (A + A');
    
    % linear term
    b = -2 * x;
    
    % minimize cost function
    y = - 2 * A \ b;
end
