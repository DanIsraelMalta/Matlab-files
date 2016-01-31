%{
% 1D integration of function handle func in region [a, b] divided to n+1
% segments, utilizing chebyshev polynomial approximation.
%
% func - function handle of one argument
% [a, b] - integration region
% n - number of segments in region [a, b]
% isVector - 'true' if f is vecorized, 'false' - otherwise
%
% example:
 
% functio to be integrated (vectrized)
F = @(x)1./(x.^3-2*x-5);
 
% matlab built in quad
tic; Q = quad(F,0,2); toc;
 
% chebyshev quad
tic; Q1 = chebyQuad(F, 0, 2, 14, 1); toc;
 
disp(['relative error to matlab quad: ', num2str(1-Q/Q1), '%']);;
 
% Dan I. Malta 2016
%}
function I = chebyQuad(func, a, b, n, isVector)
    % housekeeping
    n = n-1;
    x = cos(pi * (0:n)' / n);   % chebyshev points
 
    % is func vectorised?
    if isVector
        fx = feval(func, 0.5 * (b - a) * x + 0.5 * (b + a));
        fx = reshape(fx, numel(x), 1);
    else
        fx = zeros(size(x));
        for i = 1 : numel(fx)
            fx(i) = feval(func, 0.5 * (b - a) * x(i) + 0.5 * (b + a));
        end
    end
 
    % Chebyshev coefficients
    fx = fx / (2 * n);
    g = real(fft(fx([1 : (n + 1), n : (-1) : 2])));
    c = [g(1); g(2 : n) + g(2 * n : (-1) : (n + 2)); g(n + 1)];
    w = zeros(size(c'));
    w(1 : 2 : end) = 2 ./ (1 - (0 : 2 : n).^2);
 
    % output
    I = 0.5 * (b - a) * (w * c);
end
