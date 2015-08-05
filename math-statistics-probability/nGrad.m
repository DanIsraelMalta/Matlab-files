%{
% nGrad calculates, numerically, the gradient vector of a given vectorial function
%
% fcn - function (accpts nx1 vector input)
% x0 - gradient location (nx1 vector)
% dx - gradient step size
% type - gradient type:
%        'forward' - calculate the gradient of f(x0+dx) - f(x0)
%        'centeral' - calculate the gradient of f(x0+dx) - f(x0 - dx)
% J - gradient vector
%
% Dan I. Malta 2008
%}
function g = nGrad(fcn, x0, dx, type)
  % housekeeping
  dx = ones(size(x0)) * dx;
  N = length(x0);

  % gradient
  if strcmp(type, 'forward')
    f0 = fcn(x0);
    M = length(f0);
    g = zeros(N, 1);
    for j = 1 : N
      xj = x0;
      xj(j) = x0(j) + dx(j);
      fj = fcn(xj);
      g(j) = (fj - f0) / dx(j);
    end
  elseif strcmp(type, 'centeral')
    for j = 1 : N
      xj1 = x0;
      xj2 = x0;
      xj1(j) = x0(j) - dx(j);
      xj2(j) = x0(j) + dx(j);
      fj1 = fcn(xj1);
      fj2 = fcn(xj2);
      if j == 1
        M = length(fj1);
        g = zeros(N, 1);
      end
      g(j) = (fj2 - fj1) / 2 / dx(j);
    end
  end
end
