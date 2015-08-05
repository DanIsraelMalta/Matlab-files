%{
% levinsonDurbin - performs the levinson-durbin recrusion
%
% R - vector holding autocorrelation lags [R(0), R(1), ... , R(M)] 
% L - unit lower triangular matrix holding the reversed prediction error filters
% E - mean squared prediction error vector [E(0), E(1), ... , E(M)] 
%
%
% Dan I. Malta 2010
%}
function [L, E] = levinsonDurbin(R)
    % housekeeping
    M = length(R) - 1;
    R = R(:);
    L = eye(M + 1);
    E(1) = R(1);
    
    % update L rows recrusively
    for p = 1 : M
        gamma = L(p, 1 : p) * R(2 : p + 1) / E(p);
        L(p + 1, 1) = -gamma;
        if p >= 2 
            L(p + 1, 2 : p) = L(p, 1 : p - 1) - gamma * conj(flip(L(p, 1 : p - 1)));  
        end
        E(p + 1) = E(p) * (1 - abs(gamma)^2);
    end
end
