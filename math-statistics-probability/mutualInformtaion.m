%{
% mutualInformtaion estimates the mutual infortmation (and mutual/marginal probability density) between two vectors x & y.
%
% remarks:
% * maintaing connection to shanon entropy, mutual information logarithm is
%    in base 2 (so we mutual informtaion returns is in bits).
% * although not optimal, bin size (for information and density
%    estimation purposes) is fixed for both x & y.
% * this estimation schema requires x & y density's to be propr (not necessarily smooth)
%
% x, y - two vectors (must be real & non constant)
% px - estimated marginal probability density of x
% py - estimated marginal probability density of y
% pxy - estimated mutual probability density of x & y
% mi - mutual information between x & y (logarithm is in base 2)
%
% example:

x = rand(1,3000);
x1 = 0.5 + x;
x2 = 0.75 + x;
x3 = 2 + x;

[px, py, pxy, z0] = mutualInformtaion(x, x);
[px, py, pxy, z1]  = mutualInformtaion(x, x1);
[px, py, pxy, z2]  = mutualInformtaion(x, x2);
[px, py, pxy, z3]  = mutualInformtaion(x, x3);

%
% Dan I. Malta 2014
%}
function [px, py, pxy, mi] = mutualInformtaion(x, y)
    % housekeeping
    n = length(x);
    
    % bin'ing (space partition)
    numOfBins = max(floor(sqrt(n / 2)), 100);
    xyMin = min(min(x), min(y));
    xyMax = max(max(x), max(y));
    xyStep = (xyMax - xyMin) / numOfBins;
    edge = xyMin : xyStep : xyMax;
    
    % preallocation
    nx = zeros(numOfBins + 1, 1);
    px = nx;
    ny = nx';
    py = ny;
    pxy = zeros(numOfBins + 1, 1);
    
    % columnize x and rowize y
    x = x(:);
    y = reshape(y, 1, n);

    % binned marginal probability density
    nx = histc(x, edge);
    ny = histc(y, edge);
    px = nx(:) / n;
    py = ny(:) / n;
    
    % binned mutual probability density
    if min(x) < max(y) || min(y) < max(x)
        mut = px & py;      
        pxy(mut) = min(px(mut), py(mut));
    end

    % mutual information estimation
    mi = real(pxy .* log2(pxy ./ (px .* py)));
    mi(isnan(mi) | isinf(mi)) = [];
    mi = sum(mi);
end
