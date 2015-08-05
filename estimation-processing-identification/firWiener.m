%{
% firWiener  design a finitie impulse reponse wiener filter using
% levinson-durbin algorithm
%
% R - autocorrelation lag vector of y (Ryy(0), ..., Ryy(m))
% r - cross correlation lags between x & y (Rxy(0), ..., Rxy(m))
% h - wiener filter direct form coefficients
% g - wiener filter lattice coefficients
% gamma - lattice reflection coefficients (dimension is smaller by 1 from h & g vectors)
%
% Dan I. Malta 2011
%}
function [h, g, gamma] = firWiener(R, r)
    % columnize
    r = r(:);
    
    % levinson-durbin
    [L, E] = levinsonDurbin(R);
    
    % L first column
    gamma = -L(2 : end, 1);
    
    % lattice filter coefficients
    g = diag(E) \ L * r;
    
    % direct form coefficients
    h = L' * g;
end
