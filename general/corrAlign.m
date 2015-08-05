%{
% corrAlign accepts two vectors and return the second vector after it was
% circularly shifted such that the inner product between itself and the first
% vector is maximized.
%
% x - input vectors
% y - vector to be alligned with x
% z - y vector after circularly shifted
% l - amount of circular shift, in number of cells
%
% example:

t = 0 : 0.1 : 2*pi;
x = sin(t);
y = cos(t);
[z, l] = corrAlign(x, y);
figure;
plot(t, x, 'k', t, y, 'b', t, z, '-r');
legend('x = sin', 'y = cos', ['y shifted ', num2str(l), ' samples to be alligned with x']);
grid on;

%
% Dan I. Malta 2013
%}
function [z, l] = corrAlign(x, y)
    % columnize
    x = x(:);
    y = y(:);
           
    % align
    c = ifft(conj( fft(x)) .* fft(y));
    [m, l] = max(abs(c)); %#ok<ASGLU>
    z = circshift(y, l);
end
