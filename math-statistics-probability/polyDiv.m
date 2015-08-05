%{
% given a(x) & b(x), calculate q(x) & r(x) in the following relation: a(x) = q(x) * b(x) + r(x)
%
% Dan I. Malta 2012
%}
function [q, r] = polyDiv(a, b)
    % housekeeping
    na = length(a);
    nb = length(b);
    
    % factorization
    if na >= nb && nb > 1
        nq = na - nb + 1;
        q = zeros(1, nq);
        for i = 1 : nq
            q(i) = a(i) / b(1);
            a(i : i + nb - 1) = a(i : i + nb - 1)  - q(i) * b;
        end
    elseif nb == 1
        q = a / b;
        r = 0;
    else
        q = 0;
        r = a;
    end
end
