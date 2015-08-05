%{
% bootstrap estimates the percentile interval of any statistical calculated function,
% for a given confidence interval..
% bootstrap uses the bootstrap estimation method.
%
% x - data vector
% func - statistical function given as an inline (string) on variable 'x'
% alpha - confidence interval (1 - alpha) * 100% (alpha > 0.5)
% b - number of bootstrap replicates (shold be more then 1000)
% v - func output
% lo - percentile low limit
% hi - percentile high limit
%
% example:

clc;
alpha = 0.05;
x = randn(1,50);
[value, lo, hi] = bootstrap(x, 'mean(x)', alpha, 2000);
disp(['for ', num2str(length(x)), ' samples out of N(0,1), the mean value is ', num2str(value), ' (for confidence of ', num2str((1-alpha)*100), ' this calculation boundary is [', num2str(lo), ', ', num2str(hi),']']);
alpha = 0.05;
x = randn(1,1000);
[value, lo, hi] = bootstrap(x, 'mean(x)', alpha, 2000);
disp(['for ', num2str(length(x)), ' samples out of N(0,1), the mean value is ', num2str(value), ' (for confidence of ', num2str((1-alpha)*100), ' this calculation boundary is [', num2str(lo), ', ', num2str(hi),']']);
alpha = 0.05;
x = randn(1,3000);
[value, lo, hi] = bootstrap(x, 'mean(x)', alpha, 2000);
disp(['for ', num2str(length(x)), ' samples out of N(0,1), the mean value is ', num2str(value), ' (for confidence of ', num2str((1-alpha)*100), ' this calculation boundary is [', num2str(lo), ', ', num2str(hi),']']);

%
% Dan I. Malta 2015
%}
function [v, lo, hi] = bootstrap(x, func, alpha, b)
    % housekeeping
    x = x(:);
    n = length(x);
    f = inline(func);
    boot = zeros(1, b);
    
    % statistics observed value
    v = f(x);
    
    % resampling & bootstrap replicates (uniform)
    for i = 1 : b
        ind = ceil(n .* rand(n, 1));
        xstar = x(ind);
        boot(i) = f(xstar);
    end
    
    % percentile interval
    k = floor((b + 1) * alpha / 2);
    s = sort(boot);
    lo = s(k);
    hi = s(b + 1 - k);
end
