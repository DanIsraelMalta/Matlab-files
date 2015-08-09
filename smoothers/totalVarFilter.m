%{
% totalVarFilter performs the total variance filtering schema
% this smoother is especially usefull for piecewise linear / polynominal data corrupted by noise.
%
% x - signal
% lambda - smoothign parameter (larger lambda - larger smooth)
% n - number of iterations
% y - output
%
% example:

t = 0 : 1 : 200;
y(1:50) = 2 * t(1:50);
y(50:120) = -0.003 * t(50 : 120).^2 + 1.1 * t(50 : 120);
y(120 : 201) = - 2 .* t(t(120 : 201));
y = y + 20 * randn(size(t));
z = totalVarFilter(y, 200, 200);

figure;
plot(t, y, 'r', t, z, 'b');
legend('signal', 'total variance filtered signal');
grid on;

% Dan I. Malta 2011
%}
function y = totalVarFilter(x, lambda, n)
    % housekeeping
    len = length(x);
    z = zeros(1, len - 1);
    
    % objective parameters
    alpha = 4;
    T = lambda/2;
    
    % iteratively minimize J(k) = sum(abs(y-x).^2) + lambda*sum(abs(diff(y)));
    % to find the minimal total variance approximation of x
    for k  =  1 : n
        % y-D'z
        y = x - [-z(1), -diff(z), z(end)];
        
        % z(n) = z(n-1) + Dz/alpha
        z = z + 1 / alpha * diff(y);
        
        % clip(z,T)
        z = max(min(z, T), -T);
    end
end
