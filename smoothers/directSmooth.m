%{
% Direct finith difference N'th order smoother
% 
% x         = observation
% beta = smoothing parameter (the higher beta - the smoother the output)
% n         = linear smoother  (difference) order (values should be 1, 2 or 3)
% z         = smoothed observaion (minimizing the sum: (x[i]-z[i])^2 + beta * (z[i]-z[i-1])^2)
% 
% example:

t = 0 : 0.01 : 2*pi;
y = 3 * sin(t) + 2.5 * cos(t);
x = y + randn(size(t));
z = directSmooth(x, 100, 1);

plot(t,x, 'r', t, y, 'b', t, z, 'k');
 grid on;
legend('observation', 'signal', 'smooth');

%
% Dan I. Malta 2014
%}
function z = directSmooth(x, beta, n)
    len = length(x);
    P = diff(eye(len), n);
    z = (speye(len) + beta * P' * P)  \ x';
end
