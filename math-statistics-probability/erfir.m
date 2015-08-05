%{
% erfir calculates the error function of pure imaginary or real argument
%
% z - argument
% f - error function of z
%
% Dan I. Malta 2008
%}
function f = erfir(z)
    % housekeeping
    y = imag(z);
    
    % error function
    if abs(y) < eps % "real"
        f = erf(real(z));
    elseif abs(y) <= 8   % small number approximation
        n = 1 : 32;
        f = (i / pi) * (2 * sum(exp(-0.25 * n .* n) ./ n .* sinh(n * y))+y); 
    elseif abs(y) < 27 % medium approximation
        m = 193;
        s = 1;
        y2 = 2 * y^2;
        for n = m : -2 : 1
            s = 1 + n * s / y2;
        end
        f = i * s * exp(y^2) / (y * sqrt(pi));
    else
        f = sign(y) * i * inf;
    end
end
