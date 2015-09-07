%{
% danSgolay performs the savitsky-golay smoothing mehtod upon a uniformely spaced data set.
%
% x, y - data set to be smoothed
% win - span of smoothing window [samples]
% ord - order / degree of polynomial fit
% ys - smoothed y
%
% example:

t = 0 : 0.01 : 2*pi;
y = 3 * sin(t) + 2.5 * cos(t);
x = y + randn(size(t));
z = danSgolay(t, x, 50, 2);

plot(t,x, 'r', t, y, 'b', t, z, 'k');
 grid on;
legend('observation', 'signal', 'smooth');

%
% Dan I. Malta 2015
%}
function ys = danSgolay(x, y, win, ord)
    % columnize
    x = x(:);
    y = y(:);
    
    % housekeeping
    n = length(x);
    win = win + rem(win - 1, 2);  % window must be odd
    halfWin = (win - 1) / 2;
    
    % polynomial regression base
    v = ones(win, ord + 1);
    xw = (-halfWin : halfWin)';
    for i = 1 : ord
        v(:, i + 1) = halfWin .^ i;
    end
    [q, r] = qr(v, 0);
    
    % savitsky-golay utilizing linear filter
    y2 = filter(q * q(halfWin + 1, :)', 1, y);
    
    % edge handling (...should be improved...)
    y1 = q(1 : halfWin, :) * q' * y(1 : win);
    y3 = q((halfWin + 2) : end, :) * q' * y(n - win + 1 : n);
    
    % output
    ys = [y1; y2(win : end); y3];
end
