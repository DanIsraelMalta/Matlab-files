%{
% logProb calculates the log(P(lo < Z < hi)), where Z~N(0,1), for the region [lo, hi]
%
% Remarks:
% - calculation is percise and apathetic for rounding errors, but it is slow.
% - region borders order is not importance ([lo, hi], [hi, lo]), logProb takes care of it.
%
% lo, hi - region borders
% P - log(P(lo < Z < hi))
%
% Dan I. Malta 2015
%}
function P = logProb(lo, hi)
    % log pf tail of Z~N(0, 1)
    logPhi = inline('-0.5 * x.^2 - log(2) + reallog(erfcx(x / sqrt(2)))');
    
    % probability logarithm calculation
    if lo > 0                                   % hi > lo > 0
        Plo = logPhi(lo);    % upper tail log
        Phi = logPhi(hi);    % lower tail log
        P = Plo + log1p(-exp(Phi - Plo));
    elseif  hi < 0                         % lo < hi < 0
        Plo = logPhi(-lo);    % upper tail log
        Phi = logPhi(-hi);    % lower tail log
        P = Plo + log1p(-exp(Plo - Phi));
    elseif (lo < 0) && (hi > 0)     % lo < 0 < hi
        Plo = erfc(-lo / sqrt(2)) / 2; % lower tail
        Phi = erfc(hi / sqrt(2)) / 2;  % upper tail
        P = log1p(-Plo - Phi);
    else
        P = 0;
    end
    P = real(P);
end
