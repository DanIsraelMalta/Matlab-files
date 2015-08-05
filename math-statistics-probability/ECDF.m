%{
% ECDF calculates a vector empirical cumulative distribution function
%
% x - vector of observation
% cdf - empirical CDF
%
% example:

x = randn(1, 3000) / 3;
y = 2 * rand(1, 1000) - 1;
cdfx = ECDF(x);
cdfy = ECDF(y);

figure;
plot(sort(x), cdfx, 'b', sort(y), cdfy, 'r');
legend('normal cdf', 'uniform cdf');
xlabel('x');
ylabel('CDF');
grid on;

%
% Dan I. Malta 2008
%}
function cdf = ECDF(x)
    bins = sort(x);
    counts = histc(x, bins, 1);
    cdf = cumsum(counts) / sum(counts);
end
