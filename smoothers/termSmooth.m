%{
% A specially designed smoother to suppress derivative  by adjusting the
% impulse transfer function
%
% x         = obseravation
% tau    = smoothign factor (the higher tau -> the smoother the outcome)
% type = smoothing type:
%                0 - output has a non negative impulse reponse
%                0 - output has a perfect response to gradual curvature
% y         = smoothed output
% example:

t = 0 : 0.01 : 2*pi;
y = 3 * sin(t) + 2.5 * cos(t);
x = y + randn(size(t));
z0 = termSmooth(x, 70, 0);
z1 = termSmooth(x, 70, 1);

plot(t,x, 'r', t, y, 'g', t, z0, 'k', t, z1, 'b');
 grid on;
legend('observation', 'signal', 'smooth #0 (non negative impulse reponse)', 'smooth #1 (perfect response to gradual curvature)');

%
% Dan I. Malta (2014)
%}
function y = termSmooth(x, tau, type)
switch type
    case 0
        n = 1 + 2 * floor(1.990 * tau / 2);                  % odd filter length
        t = (pi / (n+1)) * ((1-n) : 2 : (n-1))';                    % special linspace
        h = 0.8144 + cos(t) + 0.1856 * cos(2*t);          % 3-term filter
    case 1
        n = 1 + 2 * floor(3.741 * tau / 2);                  % odd filter length
        t = (pi / (n+1)) * ((1-n) : 2 : (n-1))';                    % special linspace
        h = 0.8855 + cos(t * (1:3)) * [1.684 .9985 .2]';   % filter with flat-top response
end

% filter
y = conv(x, h) ./ conv(ones(size(x)), h);
y = y((n + 1) / 2 : end - (n - 1) / 2);
end
