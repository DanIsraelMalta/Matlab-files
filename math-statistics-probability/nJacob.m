%{
% nJacob calculates, numerically, the jacobian matrix of a given vectorial function
%
% fcn - function (accpts nx1 vector input)
% x0 - jacobian location (nx1 vector)
% dx - jacobian step size
% type - jacobian type:
%        'forward' - calculate the jacobian of f(x0+dx) - f(x0)
%        'centeral' - calculate the jacobian of f(x0+dx) - f(x0 - dx)
% J - jacobian matrix
%
% Dan I. Malta 2008
%}
function J = nJacob(fcn, x0, dx, type)
  % housekeeping
  dx = ones(size(x0)) * dx;
  N = length(x0);
  
  % jacobian
  if strcmp(type, 'forward')
    f0 = fcn(x0);
    M = length(f0);
    J = zeros(M, N);
    for j = 1 : N
      xj = x0;
      xj(j) = x0(j) + dx(j);
      fj = fcn(xj);
      J(:, j) = (fj - f0) / dx(j);
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
        J = zeros(M, N);
      end
      J(:, j) = (fj2 - fj1) / 2 / dx(j);
    end
  end
end
