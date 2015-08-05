%{
% given the polynomial fraction
 %    c(x)
% --------------     a(x) degree is p, b(x) degree is q, c(x) degree is p+q-1
% a(x) * b(x)
%
% calculate the following partial fraction:
% s(x)      t(x)
% ------ + -------
% a(x)       b(x)
%
% Dan I. Malta 2012
%}
function [s, t] = partialFrac(a, b, c)
    % housekeeping
    p = length(a) - 1;
    q = length(b) - 1;
    A = zeros(p + q);
    
    % matrix form
    for i = 1 : p + q
        if i <= p
            A(:, i) = [zeros(i - 1, 1); b'; zeros(p - i, 1)];
        else
            j = i - p;
            A(:, i) = [zeros(j - 1, 1); a'; zeros(q - j, 1)];
        end
    end
    
    % c factorization
    if length(c) < p + q
        c = [zeros(1, p + q - length(c), c)]; %#ok<NBRAK>
    end
    
    % factorization
    st = inv(A) * c';
    s = [0, st(1 : p)'];
    t = st(p + 1 : p + q)';
end
