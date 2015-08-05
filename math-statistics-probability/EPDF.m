%{
% EPDF estimate a given vector probability density using an averaged
% shifted histogram method.
%
% x - vector
% h - window width
% m - number of shifted histograms
% pd - estimated probability density
% bin - bins for fhat
%
% example:

x = 3 + randn(1, 9000);
[xpd, xbin] = EPDF(x, 1, 300);

figure;
stairs(xbin, xpd);
xlim([0 6]);
title('N(3, 1) probability density')
grid on;

%
% Dan I. Malta 2010
%}
function [pd, bin] = EPDF(x, h, m)
    % housekeeping
    n = length(x);
    delta = h / m;

    % Get the mesh.
    t0 = 0;
    tf = 20 + max(x);
    nbin = ceil((tf - t0) / delta);
    bin = t0 : delta : (t0 + delta * nbin);
    
    % bin counts for smaller binwidth delta
    vk = histc(x, bin);
    pd = [zeros(1, m - 1), vk, zeros(1,m - 1)];
    
    % weight vector.
    kern = inline('(15 / 16) * (1 - x.^2).^2');
    ind = (1 - m) : (m - 1);
    den = sum(kern(ind / m));
    wm = m * (kern(ind / m)) / den;
    
    % bin heights over smalle bins
    pdk = zeros(1, nbin);
    for k = 1 : nbin
        ind = k : (2 * m + k - 2);
         pdk(k) = sum(wm .* pd(ind));
    end
    pdk = pdk / (n * h);
    bc = t0 + (( 1 :k) - 0.5 ) *delta;

    pd = [pdk pdk(end)];
end
