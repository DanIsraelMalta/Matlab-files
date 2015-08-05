%{
% splineSmooth returns the value of the smoothing spline of a given data series.
%
% why not use matlab builtin spline function? well, because...
% 1) matlab spline doesn't allow you to choose the spline smoothing parameter, it automatically choose it for you.
% 2) in order to find the spline coefficients, matlab spline function solves a tridiagonal linear system
%      while I use the cross-validation method with sparse matrix. Although it is slower, it is less sensitive to numerical errors...
% 3) matlab spline performs interpolation, while I just want a smoother.
%
% x, y       - data series
% alpha - spline smoothing parameter (larger alpha -> smoother spline)
% xo         - unique values of x
% yhat    - spline values at x0
% S            - spline smoothing matrix. usefull for cross-validation purposes
%
% example:

x = 1 : 0.5 : 52;
y = 8 + 0.8 * sin(x + rand()) - 0.002* x.^2 + randn(size(x));
[x0,yhat,s] = splineSmooth(x,y,5);
yy = spline(x,y,x);
plot(x,y,'.b',x0,yhat, 'k', x, yy, 'r')
legend('obseration', 'splineSmooth', 'matlab builtin spline');
grid on;

%
% Dan I. Malta 2010
%} 
function [x0,yhat,S] = splineSmooth(x,y,alpha)
    % housekeeping
    if alpha <= 0
        error('Alpha must be greater than 0.') %#ok<ERTAG>
    end

    % sort x values.
    [xs,ind] = sort(x);
    ys = y(ind);
    
    % columnize
    x = xs(:); 
    y = ys(:);
    n = length(x);
    
    % weight vector (initialized to 1)
    w = ones(n,1);
    
    % remove "tied observation" and replace them with their average value
    h = diff(x);
    ind0 = find(h==0);
    if ~isempty(ind0)
        xt = x;
        yt = y;
        wt = w;
        i = 1;
         while ~isempty(ind0)
             indt = find(x(ind0(end)) == x);
             ym(i) = mean(y(indt));
             xm(i) = x(indt(1));
             wm(i) = length(indt);
             i = i+1;
             xt(indt) = [];
             yt(indt) = [];
             wt(indt) = [];
             [c, ia, ib] = intersect(indt, ind0); %#ok<ASGLU>
             ind0(ib) = [];
         end
         
          xu = [xt(:); xm(:)];
          yu = [yt(:); ym(:)];
          wu = [wt(:); wm(:)];
          [xus, inds] = sort(xu);
          yus = yu(inds);
          wus = wu(inds);
          x = xus;
          y = yus;
          w = wus;
    end
     n = length(x);
     
     % h[i] = x[i+1] - x[i]
     h = diff(x);
     
     % 1/h[i]
     hinv = 1 ./ h;
     W = diag(w);
     
     % my way of building the spline sparse matrix (slow but numerically stble)
     qDs = -hinv(1 : n - 2) - hinv(2 : n - 1);
     I = [1 : n - 2, 2 : n - 1, 3 : n];
     J = [2 : n - 1, 2 : n - 1,2 : n - 1];
     S = [hinv(1 : n - 2), qDs, hinv(2 : n - 1)];
     
     % Create a sparse matrix
     Q = sparse(I, J, S, n, n);
     
     % Delete the first and last columns
     Q(:, n) = [];
     Q(:, 1) = [];
     
     % find the R matrix. 
     I = 2 : n - 2;
     J = I + 1;
     tmp = sparse(I, J, h(I), n, n);
     t = (h(1 : n - 2) + h(2 : n - 1)) / 3;
     R = tmp' + tmp + sparse(2 : n - 1, 2 : n - 1, t, n, n);
     
     % Get rid unwanted rows/columns
     R(n, :) = [];
     R(1, :) = [];
     R(:, n) = [];
     R(:, 1) = [];
     
     % spline
     S1 = Q' * y;
     S2 = R + alpha * Q' * inv(W) * Q;
     
     % Solve for gamma
     gam = S2 \ S1;
     
     % fhat
     yhat = y - alpha * inv(W) * Q * gam;
     S = inv(W + alpha * Q * inv(R) * Q') * W;
     x0 =x; 
end
