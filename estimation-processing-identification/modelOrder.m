%{
% modelOrder - estimates a given signal subspace order (for future estimation purposes) using AIC & MDL criteria's
%
% N - length of observation vector
% lambda - vector holding signal sample covariance eigenvalues, in increasing order
% aic - vector holding AIC[k]
% mdl - vector holding MDL[k]
%
% notes:
% - choose minimum value of aic/mdl for optimal estimated model order
% - although MDL is asymptotically consistnen, AIC is not...
%
% Dan I. Malta 2010
%}
function [aic,mdl] = modelOrder(N,lambda)
    % housekeeping
    M = length(lambda) - 1;
    
    % criteria calculation
    for k = 1 : M + 1
        la = lambda(1 : k);
        L = mean(log(la)) - log(mean(la));
        aic(k) = -2 * N * k * L + 2 * (M + 1 - k) * (M + 1 + k);
        mdl(k) = -N * k * L + 0.5 * (M + 1 - k) * (M + 1 + k) *log(N);
    end
