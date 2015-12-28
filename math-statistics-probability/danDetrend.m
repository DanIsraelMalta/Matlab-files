%{
% danDetrend removes the trend from a given vector
%
% why use danDetrend and not matlab builtin detrend function? because:
% 1) danDetrend handles NaN's
% 2) danDetrend can handle unequally spaced data
% 3) danDetrend returns, along with the detrended data, the removed trend
% 4) danDetrend can detrend not only constant and linear trends, but also
%     polynoomial trend of any order.
%
% [...] = danDetrend(x, 0) -> removes constant mean
% [...] = danDetrend(x, 1) -> removes a linear trend
% [...] = danDetrend(x, p) -> removes a polynomial trend of order p
% [...] = danDetrend(x, t, p) -> removes a polynomial trend of order p from  series {t, x}
%
% [T, X]=danDetrend(...)  -> X = the detrended data,  T = the removed trend
%
% example:

% parabula detrend!!!
x = -10 : 0.1 : 10;
y = x.^2 + 7*randn(size(x));
[trend, y_detrend] = danDetrend(y, 2);

figure;
subplot(2,1,1);
plot(x, y, 'b', x, y_detrend, 'r');
legend('white noise with parabula trend', 'estimated trend');
grid on;
subplot(2,1,2);
plot(x, trend);
title('white noise with after parabula trend was removed -> can be analyzed as wide-sense stationary!');
grid on;

%
% Dan I. Malta 2011
%}
function [X, T] = danDetrend(t, X, p)
    % housekeeping
    if (nargin == 1)
        p = 1;
        X = t;
        t = [];
    elseif (nargin == 2)
        if strcmpi(X, 'constant')
            p = 0;
            X = t;
            t = [];
        elseif strcmpi(X, 'linear')
            p = 1; 
            X = t; 
            t = [];
        elseif all(size(X) == 1)
            p = X;
            X = t;
            t = [];
        else
            p = 1;
        end
    end

    % columnize data
    [m, n] = size(X);
    if (m == 1)
        X = X';
        r = n;
    else
        r = m;
    end

    % time scale check
    if isempty(t),
        t = (1 : r).'; % create time scale
    elseif ~all(size(t) == size(X)) 
        t = t(:);
    end
    
    % check t & X dimesions
    if ~all(size(X, 1) == size(t, 1))
        fprintf (2,'danDetrend: size(t,1) must same as size(x,1) \n');
    end
    
    % check the order of the polynomial 
    if (~(all(size(p)==1) & (p == round (p)) & (p >= 0)))
        fprintf (2,'danDetrend:  p must be a nonnegative integer\n');
    end

    % do we need more memory?
    if nargout > 1
        T = zeros(size(X)) + nan;       
        
        % for multiple time scales
        if size(t, 2) > 1
            for k = 1 : size(X,2)
                idx = find(~isnan(X(:, k)));
                b = (t(idx, k) * ones(1, p + 1)) .^ (ones(length(idx), 1) * (0 : p));
		        T(idx, k) = b * (b \ X(idx, k));
            end
        else			% if only one time scale is used
            b = (t * ones(1, p + 1)) .^ (ones(length(t),1) * (0 : p));
            for k = 1 : size(X, 2)
                idx = find(~isnan(X(:, k)));
                T(idx,k) = b(idx, :) * (b(idx, :) \ X(idx, k));
            end
        end
        
        % detrend
        X = X - T;
        
        if (m == 1)
	        X = X';
	        T = T';
        end
    else % needs less memory
        % for multiple time scales
        if size(t,2) > 1
            for k = 1 : size(X,2)
                idx = find(~isnan(X(:, k)));
                b = (t(idx, k) * ones(1, p + 1)) .^ (ones(length(idx), 1) * (0 : p));
                X(idx, k) = X(idx, k) -  b * (b \ X(idx, k));
            end
        else			% if only one time scale is used
            b = (t * ones(1, p + 1)) .^ (ones(length(t),1) * (0 : p));
            for k = 1 : size(X, 2)
                idx = find(~isnan(X(:, k)));
                X(idx, k) = X(idx, k) - b(idx, :) * (b(idx, :) \ X(idx, k));
            end
        end
        
        if (m == 1)
	        X = X';
        end
    end
end
