%{
% confCorr calculates the confidence intervals of th correlation coefficient
%
% requires statistics toolbox.
%
% r - correlation
% n - sample size (must be larger then 2)
% alpha - desired confidence level (0.05 will yoeld a 1-alpha = 0.95 = 95% confidence)
% cl - lower and upper confidence bounds
% p - p value of correlation bi-directional test
%
% example:

% how good can we estimata a correlation value of 0.5 calculated from sample size of 20?
[ci, p] = confCorr(0.5, 20, 0.05);

%
% Dan I. Malta 2014
%}
function [cl, p] = confCorr(r, n, alpha)
    % correlation to Fisher's z
    fz = log((1 + r)/(1 - r)) / 2;
    
    % standard errpr of fisher's's Z-score
    ser = 1 / sqrt(n - 3);
    Zcrit = norminv(alpha / 2, 0, 1);
    
    % confidence bound
    lolim = fz - abs(Zcrit) * ser;
    uplim = fz + abs(Zcrit) * ser;
    lolim_r = (exp(1) ^ (2 * lolim) - 1) / (exp(1) ^ (2 * lolim) + 1);
    uplim_r = (exp(1) ^ (2 * uplim) - 1) / (exp(1) ^ (2 * uplim) + 1);
    cl = [lolim_r uplim_r];
    
    % non directional p-value
    t = r * sqrt((n - 2) / (1 - r^2));
    p = 2 * (1 - tcdf(t, n - 2));
end
