%{
% adaptWiener implements an adaptive wiener filter structure
%
%  given a pair of input samples {x, y}, calculate the pair {xhat, e} and
%  update internal state w, and adapt filter weights h.
%
% prior to calling adaptWiener, initialize as follows:
%   w = zeros(M + 1, 1);
%   h = zeros(M + 1, 1);
%   a = [1; zeros(M, 1)];
%   b = [zeros(M, 1); 1];
%   D1a = 0.01;
%   D1b = 0.01;
%   kt = zeros(M, 1);
%   nt = 0;
%
% x - input sample (primary)
% y - input sample (secondary)
% w - internal stat vector (w(0), ..., w(m))
% h - wiener filter coefficients vector (h(0), ..., h(m))
% a - apriori forward prediction filter
% b - apriori backward prediction filter
% Dla - forward prediction error (from previos sample)
% Dlb - backward prediction error (from previos sample)
% kt - kalman gain vector
% nt - likelihood variable
% lambda - "forgetting" factor (fadin memeory factor)
%
% xhat - estimated x
% w - updated delay-linear vector
% h - updated wiener filter coefficients
% a - updated aposteriori forward prediction filter
% b - updated aposteriori backward prediction filter
% Dla - updated forward prediction error (from previos sample)
% Dlb - updated backward prediction error (from previos sample)
% kt - updated kalman gain vector
% nt - updated likelihood variable
%
% Dan I. Malta 2012
%}
function [xhat, e, w, h, a, b, D1a, D1b, kt, nt] = adaptWiener(x, y, w, h, a, b, D1a, D1b, kt, nt, lambda)
    % initialization
    M = length(h) - 1;
    w(1) = y;
    e0a = a' * w;
    e1a = e0a / (1 + nt);
    D0a = lambda * D1a;
    k(1) = pinv(D0a) * e0a;
    D1a = D0a + e1a * e0a;
    k = [0; kt] + k(1) * a;
    D0b = lambda * D1b;
    e0b = k(M + 1) * D0b;
    kbar = k(1 : M) - k(M + 1) * b(1 : M);
    
    % ahead calculation
    nu = nt + e0a * k(1, 1);
    nubar = nu - e0b * k(M + 1);
    e1b = e0b / (1 + nubar);
    D1b = D0b + e1b * e0b;
    
    % initial calculation
    a = a - e1a * [0; kt];
    b = b - e1b * [kbar; 0];
    
    % wiener
    xhat0 = h' * w;
    e0 = x - xhat0;
    e = e0 / ( 1 + nu);
    xhat = x - e;
    
    % update
    h = h + e * k;
    kt = kbar;
    nt = nubar;
    w(M + 1 : -1 : 2) = w((M + 1 : -1 : 2) - 1);
end
