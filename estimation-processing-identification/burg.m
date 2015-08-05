%{
% burg performs burg's method of linear prediction
%
% y - observation vector
% M - desired AR model order
% L - unit lower triangular matrix holding the reversed prediction error filters
% E - mean squared prediction error vector [E(0), E(1), ... , E(M)] 
%
% Dan I. Malta 2010
%}
function [L, E] = burg(y,M)
    % housekeeping
    N = length(y);
    y = y(:);
    ea = y;
    eb = y;
    L = eye(M + 1);
    E(1) = abs(y' * y) / N;                              

    % burg's
    for p = 1 : M
        n = N - 1 : -1 : p;          
        gamma = 2 * eb(n)' * ea(n+1) / (ea(n+1)' * ea(n+1) + eb(n)' * eb(n));
        
        % order-p lattice section
        temp = ea(n+1);
        ea(n+1) = temp - gamma * eb(n);
        eb(n+1) = eb(n) - conj(gamma) * temp;
        
        % L rows update
        L(p + 1, 1) = -gamma;
        if p >= 2
            L(p + 1, 2 : p) = L(p, 1 :p - 1) - gamma * conj(flip(L(p, 1 : p - 1)));  
        end
        E(p+1) = E(p) * (1 - abs(gamma)^2);
    end
end
