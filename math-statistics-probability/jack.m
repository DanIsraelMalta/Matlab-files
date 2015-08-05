%{
% jack estimates the bias and standard error of any statistical calculated function.
% jack uses the jackknife estimation method.
%
% x - data vector
% func - statistical function given as an inline (string) on variable 'x'
% v - func output
% b - func estimated bias
% s - func estimated standard error
%
% example:

clc;
x = randn(1,50);
[value, bias, std] = jack(x, 'mean(x)');
disp(['for ', num2str(length(x)), ' samples out of N(0,1), the mean value is ', num2str(value), ' (this calculation bias is ', num2str(bias), ', and standard error is ', num2str(std), ')']);
x = randn(1,3000);
[value, bias, std] = jack(x, 'mean(x)');
disp(['for ', num2str(length(x)), ' samples out of N(0,1), the mean value is ', num2str(value), ' (this calculation bias is ', num2str(bias), ', and standard error is ', num2str(std), ')']);
x = randn(1,3000);
[value, bias, std] = jack(x, 'mean(x)');
disp(['for ', num2str(length(x)), ' samples out of N(0,1), the mean value is ', num2str(value), ' (this calculation bias is ', num2str(bias), ', and standard error is ', num2str(std), ')']);

%
% Dan I. Malta 2015
%}
function [v, b, s] = jack(x, func)
    % housekeeping
    x = x(:);
    n = length(x);
    f = inline(func);
    j = zeros(size(x));
    
    % statistics observed value
    v = f(x);
    
    % jackknife
    for i = 1 : n
        j(i) = f([x(1 : i - 1); x(i + 1 : n)]);
    end
    
    % output
    b = (n - 1) * (mean(j) - v); 
    s = sqrt((n - 1) / n * sum((j - mean(j)).^2));
end
