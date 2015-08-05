%{
% loess2D performs a 2D loess regression/smoothing (alpha = 0.5, 1inear) for a given data set.
% loess2D uses loessFit.m function.
%
% x, y - spatial observation
% xo, yo - loess regression
%
% example:

t = 0 : 0.05 : 2 * pi;
x = sin(t) + 0.1 * randn(size(t));
y = cos(t) + 0.1 * randn(size(t));
[xs, ys] = loess2D(x, y);

figure;
plot(x, y, 'r', xs, ys, 'b');
legend('original', 'LOESS2D smooth');
grid on;

%
% Dan I. Malta 2009
%}
function [xo, yo] = loess2D(x, y)
    % coulmnize
    x = x(:);
    y = y(:);
    
    % housekeeping
    n = length(x);
    J = ceil(n / 2);
    MAD = inline('median(abs(x - median(x)))');
    
    % MAD normalization
    xstar = (x - median(x)) / MAD(x);
    ystar = (y - median(y)) / MAD(y);
    
    % normalized bounds
    s = ystar + xstar;
    d = ystar - xstar;
    sstar = s / MAD(s);
    dstar = d / MAD(d);
    
    % polar
    [th, m] = cart2pol(sstar, dstar);
    
    % radius
    z = m.^(2/3);
    
    % loess
    tx = -2 * pi + th((n - J + 1) : n);
    ntx = length(tx);
    tx = [tx; th];
    tx = [tx; th(1 : J)];
    ty = z((n - J + 1) : n);
    ty = [ty; z];
    ty = [ty; z(1:J)];
    tyhat = loessFit(tx, ty, tx, 0.5, 1);
    
    % cartesian
    tyhat(1 : ntx) = [];
    mhat = tyhat(1 : n).^(3/2);
    [shatstar, dhatstar] = pol2cart(th, mhat);
    
    % original scales.
    shat = shatstar * MAD(s);
    dhat = dhatstar * MAD(d);
    xhat = ((shat - dhat) / 2) * MAD(x) + median(x);
    yhat = ((shat + dhat) / 2) * MAD(y) + median(y);
    K = convhull(xhat,yhat);
    xo = xhat(K);
    yo = yhat(K);
end
