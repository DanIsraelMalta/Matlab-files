%{
% autocorrExt extendes an autocorrelation sequence using levinson durbin recrusion
%
% Rold - vector holding old autocorrelation lags
% gamma - vector holding new reflection coefficients (magnitude must be smaller then 1)
% R - vector holding extended autcorrelation lags
% L - unit lower triangular matrix holding the reversed prediction error filters
% E - mean squared prediction error vector [E(0), E(1), ... , E(M)] 
%
% notice that maximum entropy extension corresponds to gamma = [0, ..., 0];
%
% Dan I. Malta 2011
%}
function [R,L,E] = autocorrExt(Rold,gamma)
    % housekeeping
    M = length(Rold) - 1;
    q = length(gamma);
    R = Rold(:);
    
    % calculate old L & E
    [L, E] = levinsonDurbin(R);
    
    % add q rows to L
    for p = M + 1 : M + q
        L(p + 1, p + 1) = 1;                             
        L(p + 1, 1) = -gamma(p - M);                     % p-M = 1:q
        if p >= 2 
            L(p + 1, 2 : p) = L(p, 1 : p - 1) - gamma(p - M) * conj(flip(L(p, 1 : p - 1)));  
        end
        
        % next error
        E(p + 1) = E(p) * (1 - abs(gamma(p - M))^2);
        
        % next lag
        R(p + 1, 1) = - L(p + 1, 1 : p) * R(1 : p);
    end
end
