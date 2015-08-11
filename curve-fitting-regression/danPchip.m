%{
% danPchip calculates the piecewise cubic hermite interpolating polynomial,
% exactly like matlab's builtin pchip function. The only difference is that
% this function also returns first and second derivatives.
%
% x, y - 2D data series
% xi - x values for output
% type - tangent type calculation:
%               'finite' - finite difference
%               'catmull' - catmull-rom spline
%               'monotone' - monotone interpolation
% yi - interpolant at xi
% ypi - first order derivative at xi
% yppi - second order derivative at xi
%
% Dan I. Malta 2009
%}
function [yi, ypi, yppi] = danPchip(x, y, xi, type)
    % columnize
    x = x(:);
    y = y(:);
    xi = xi(:);
    
    % housekeeping
    n = length(x);
    ni = length(xi);
    yi = zeros(size(xi));
    ypi = zeros(size(xi));
    yppi = zeros(size(xi));
    
    % interpolation
    for i = 1 : ni
        % find xi location
        lo = 1;
        hi = n;
        while hi - lo > 1
            k = fix((hi + lo) / 2);
            if x(k) > xi(i)
                hi = k;
            else
                lo = k;
            end
        end
        
        % tangent calculation
        h = x(hi) - x(lo);
        if strcmp(type, 'finite')
            if lo == 1
                a = (y(hi) - y(lo)) / h;
                b = (y(hi + 1) - y(hi)) / (2 * (x(hi + 1) - x(hi))) + (y(hi) - y(lo)) / (2 * h);
            elseif hi == n
                a = (y(hi) - y(lo)) / (2 * h) + (y(lo) - y(lo - 1)) / (2 * (x(lo) - x(lo - 1)));
                b = (y(hi) - y(lo)) / h;
            else
                a = (y(hi) - y(lo)) / (2 * h) + (y(lo) - y(lo-1)) / (2 * (x(lo) - x(lo-1)));
                b = (y(hi + 1) - y(hi)) / (2 * (x(hi + 1) - x(hi))) + (y(hi) - y(lo)) / (2 * h);
            end
        elseif strcmp(type, 'catmull')
            if lo == 1
                a = (y(hi) - y(lo)) / h;
                b = (y(hi + 1) - y(lo)) / (x(hi + 1) - x(lo));
            elseif hi == n
                a = (y(hi) - y(lo-1)) / (x(hi) - x(lo - 1));
                b = (y(hi) - y(lo)) / h;
            else
                a = (y(hi) - y(lo-1)) / (x(hi) - x(lo - 1));
                b = (y(hi + 1) - y(lo)) / (x(hi + 1) - x(lo));
            end
        elseif strcmp(type, 'monotone')
            if lo == 1
                a = (y(hi) - y(lo)) / h;
                b = ((y(hi) - y(lo)) / h + (y(hi + 1) - y(hi)) / (x(hi + 1) - x(hi))) / 2;
            elseif hi == n
                a = ((y(lo) - y(lo - 1)) / (x(lo) - x(lo - 1)) + (y(hi) - y(lo)) / h) / 2;
                b = (y(hi) - y(lo)) / h;
            else
                a = ( (y(lo) - y(lo - 1)) / (x(lo) - x(lo - 1)) + (y(hi) - y(lo)) / h) / 2;
                b = ( (y(hi) - y(lo)) / h + (y(hi + 1) - y(hi))/(x(hi + 1) - x(hi))) / 2;
            end
            if ( y(hi) == y(lo) )
                a = 0;
                b = 0;
            else
                alpha = a / ( (y(hi) - y(lo)) / h);
                beta = b / ( (y(hi) - y(lo)) / h);
                if alpha == 0 || beta == 0
                    a = 0;
                    b = 0;
                else
                    if alpha * alpha + beta * beta > 9
                        tau = 3 / sqrt(alpha * alpha + beta * beta);
                        a = tau * alpha * (y(hi) - y(lo)) / h;
                        b = tau * beta * (y(hi) - y(lo)) / h;
                    end
                end
            end
        else
            error('unknown differentiation schema.'); %#ok<ERTAG>
        end
        
        % pchip
        t = (xi(i) - x(lo)) / h;
        h00 = 2*t.^3 - 3*t.^2 + 1;
        h10 = t.^3 - 2*t.^2 + t;
        h01 = -2*t.^3 + 3*t.^2;
        h11 = t.^3 - t.^2;
        yi(i) = h00 * y(lo) + h10 * h * a + h01 * y(hi) + h11 * h * b;
        
        % first derivative
        h00 = 6*t.^2 - 6*t;
        h10 = 3*t.^2 - 4*t + 1;
        h01 = -h00;
        h11 = 3*t.^2 - 2*t;
        ypi(i) = h00 * y(lo) / h + h10 * a + h01 * y(hi) / h + h11 * b;
        
        % second derivative
        h00 = 12*t - 6;
        h10 = 6*t - 4;
        h01 = -h00;
        h11 = 6*t - 2;
        yppi(i) = h00 * y(lo) / h.^2 + h10 * a / h + h01 * y(hi) / h.^2 + h11 * b / h;
    end
end
