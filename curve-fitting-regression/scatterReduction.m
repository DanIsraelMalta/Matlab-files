%{
% scatterReduction is my attempt to create sort of a scatterplot interpolant based on
% grouped first order difference smoothing.
% basically, its an attempt to reduce scattered data size while maintaing
% most of the data.
%
% x, y - scattered data to be smoothed & reduced
% n - number of bins (of which x shall be divided to)
% beta - smoothing parameter (the larger - the smoother)
% u - smoothed/reduced z, x-axis
% z - smoothed/reduced  z
%
% example:

x = 0 : 0.01: 12*pi;
y = sin(x) + 0.2 * randn(size(x));
[u, z] = scatterReduction(x, y, 40, 1);
figure;
plot(x, y, '-or', u, z, '-*b');
legend('original data set', 'smoothed/reduced data set');
grid on;

%
% Dan I. Malta 2015
%}
function [u, z] = scatterReduction(x, y, n, beta)
    % columnize
    y = y(:);
    x = x(:);
    
    % housekeeping
    c = zeros(1, n);
    s = zeros(1, n); 
    xmin = min(x);
    dx = (max(x) - xmin) / n;
    
    % grouping number and sum of observations
    for i = 1 : length(x)
        j = min(1 + floor((x(i) - xmin) / dx), n);
        c(j) = c(j) + 1;
        s(j) = s(j) + y(i);
    end
        
    % first order difference smoothing
    p = diff(speye(n), 1);
    z = (diag(c) + beta * p' * p) \ s';
    
    % x-axis
    u = xmin - dx/2 + (1:n) * dx;
end
