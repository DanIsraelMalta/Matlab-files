%{
% quant returns the sampled quantiles of a given vector.
%
% Unlike the statisticx toolbox (which we don't have here in building 24), I choose
% to calculate the precentiles without the usage of permuatation,
% but with the simple notion that the j-th order statistics estimates the
% (j-0.5)/n quantile -> and this observation makes the calculation extremely fast & efficient.
%
% x - observatino vector
% p - vector of required percentiles [0, 1]
% q - quantiles vector
%
% example:

x = randn(1, 3000);
q = quant(x, [0.25, 0.5, 0.75])

%
% Dan I. Malta 2007
%}
function q = quant(x, p)
    % housekeeping
    n = length(x);
    p = max(min(p, 1), 0);
    q = zeros(size(p));
    
    % quantile
    xs = sort(x);
    qh = ((1 : n) - 0.5) / n;
    q = interp1(qh, xs, p);
end
