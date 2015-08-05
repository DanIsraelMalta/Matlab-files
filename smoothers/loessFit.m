%{
% loessFit performs the loess smooth for univariate observation.
%
% This is my implementation of the LOESS smoothing/regression I used
% in T-HUMS code. It is a variant of the alogrithm found in Gelman&Rubin
% paper (1992).
% Notice that this implementation is very slow due to QR decomposition,
% but it reduces numerical errors by eliminating standard matrix inversion
% schema.
%
% x, y - spatial observation
% xo - domain (indices) over whivh output is calculated
% alpha - smoothing parameters [0, 1] (the higher - the smoother)
% deg - local regression degree (1 or 2)
% yo - loess regression @ xo
%
% example:

t = 0 : 0.01 : 2 * pi;
y = 3 * cos(t + 1.5) + 0.5 * randn(size(t));
[ts, ind] = sort(t);
ys = y(ind);
yo = loessFit(ts, ys, ts, 0.5, 2);

figure;
plot(t, y, 'r', t, 3 * cos(t + 1.5), 'k', ts, yo, 'b');
legend('observation', 'clean signal', 'LOESS smooth');
grid on;

%
% Dan I. Malta 2007
%}
function yo = loessFit(x, y, xo, alpha, deg)
    % housekeeping
    n = length(x);
    k = floor(alpha * n);
    yo = zeros(size(xo));
    
    % loess
    for i = 1 : length(xo)
        % find neighbour points
        dist = abs(xo(i) - x);
        [sdist, ind] = sort(dist);
        nxo = x(ind(1 : k));
        nyo = y(ind(1 : k));
        dxo = sdist(k);
        sdist(k + 1 : n) = [];
        
        % tricube weight function
        u = sdist / dxo;
        w = (1 - u.^3).^3;
        p = weightFit(nxo, nyo, w, deg);
        yo(i) = polyval(p, xo(i));
    end
end

% perform weighted least squares
function p = weightFit(x, y, w, deg)
    % columnize
    x = x(:);
    y = y(:);
    w = w(:);
    
    % housekeeping
    n = length(x);
    
    % matrix
    W = spdiags(w, 0, n, n);
    A = vander(x);
    A(:, 1 : length(x) - deg - 1) = [];
    
    % least squares
    v = A' * W * A;
    Y = A' * W * y;
    [q, r] = qr(v, 0);
    p = r \ (q' * Y);
    p = p';
end
