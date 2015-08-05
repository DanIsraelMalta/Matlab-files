%{
% polydeg find's the optimal degree for polynomial fitting,
% according to Akaike's information criteria (same as Bayes criteria, since we don't have prior PDF)
%
% x, y - vectors of {x, y} data
% n      - optimal polynog degree
%
% example:

x = linspace(0, 10, 300);
y = sin(x.^3 / 100) .^ 2 + 0.05 * randn(size(x));
n = optimalPoly(x,y);
ys = orthogonalPolyFit(x, y, n);
plot(x, y, '.', x, ys, 'k');
grid on;
title(['optimal degree of polynom is ', num2str(n)]);

%
% Dan I. Malta 2009
%}
function n = optimalPoly(x, y)
% array-wise
x = x(:);
y = y(:);
N = length(x);

% 0th AIC (initial guess)
p = mean(y);
ys = ones(N, 1) * p;
AIC = 2 + N * (log(2 * pi * sum((ys-y).^2) / N) + 1) + 4 / (N-2);

% search the optimal degree minimizing the Akaike's information criterion
% at least 3 steps are needed to take to remove "AIC noise" and ensure
% non local minimum
p = zeros(2,2);
p(1,2) = mean(x);
PL = ones(N,2);
PL(:,2) = x - p(1,2);
n = 1;
nit = 0;
while nit<3
    % orthogonal polynomial fitting
    if n>1
        p(1,n+1) = sum(x.*PL(:,n).^2)/sum(PL(:,n).^2);
        p(2,n+1) = sum(x.*PL(:,n-1).*PL(:,n))/sum(PL(:,n-1).^2);
        PL(:,n+1) = (x-p(1,n+1)).*PL(:,n)-p(2,n+1)*PL(:,n-1);
    end

    tmp = sum(repmat(y,[1,n+1]).*PL)./sum(PL.^2);
    ys = sum(PL.*repmat(tmp,[N,1]),2);

    % AIC (including small sample sizes correction term)
    aic = 2 * (n + 1) + N * (log(2 * pi * sum((ys-y(:)).^2) / N) + 1) + 2 * (n + 1) * (n + 2) / (N - n - 2);

    if aic>=AIC
        nit = nit+1;
    else
        nit = 0;
        AIC = aic;
    end
    n = n+1;

    if n>=N
        break
    end
end

n = n - nit - 1;
end
