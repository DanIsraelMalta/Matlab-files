%{
% qqplot draws a quantile-quantile plot (theoretical quantiles are plotted against sample order statistics)
%
% x, y - data series to be investigated
%
% example:

figure;
subplot(2,1,1);
x = randn(1, 1000);
prob = ((1:100) - 0.5) / 100;
y = -sqrt(2) .* erfcinv(2 *  prob);
qqplot(x, y);
title('normal distribution qq-plot');
subplot(2,1,2);
x = rand(1, 1000);
prob = ((1:100) - 0.5) / 100;
y = -sqrt(2) .* erfcinv(2 *  prob);
qqplot(x, y);
title('uniform distribution qq-plot');

%
% Dan I. Malta 2008
%}
function qqplot(x,y)
    % houskeeping
     m =length(x);
     n =length(y);
     xs = sort(x);
     ys = sort(y);
     
     % equal - just plot
     if m == n
         plot(xs,ys,'o')
         xlabel('X');
         ylabel('Y');
     elseif m<n % interpolate from the larger vector (y)
         % x's that go with the y's data set
         prob = ((1 : n) - 0.5) / n;
         
         % interpolation
         qs = ((1 : m) - 0.5) / m;
         ysi = interp1(prob, ys, qs, 'linear');
         
         % plot
         plot(xs,ysi,'o');
         xlabel('X');
         ylabel('Y - Interpolated');
     else % interpolate from the larger vector (x)
           %y's that go with the x's data set
         prob = ((1 : m) - 0.5) / m;
         
         % interpolation
         qs = ((1 : n) - 0.5) / n;
         xsi = interp1(prob, xs, qs, 'linear');
         
         % plot
         plot(xsi,ys,'o')
         xlabel('X - Interpolated');
         ylabel('Y');
     end
     
     grid on;
end   
