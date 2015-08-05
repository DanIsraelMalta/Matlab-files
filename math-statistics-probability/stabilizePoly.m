%{
% given a polynomial (by its coefficients), reflect all its un stable roots
% (those roots which are located outside of the unit circle) inside th
% unit circle, thuse stablizing it.
%
% a - polynomial coefficients
% b - polynomial coefficients after reflection
%
% Dan I. Malta 2009
%}
function b = stabilizePoly(a)
    % find polynomial roots (ignore zeros)
    v = roots(a);
    i = find(v ~= 0);
    
    % reflect roots in unstable
    vSign = (sign(abs(v(i)) - 1) + 1) / 2;
    v(i) = (1 - vSign) .* v(i) + vSign ./ conj(v(i));
    ind = find(a ~= 0);
    b = a(ind(1)) * poly(v);

    % make sure that real coefficients are returned real
    if ~any(imag(a))
        b = real(b);
    end
end
