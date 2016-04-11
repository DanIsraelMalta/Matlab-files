%{
% quadMin solves a constrained quadratic eigenvalue problem.
% Given a matrix A and a vector b, quadMin will minimize the L2 norm of
% |A*x-b|^2 such that x = xnorm.
%
% quadMin performs the same as matlab builtin LSQR function, but it is
% faster on certain problems which arise in the field of control &
% estimation.
%
% A - matrix (system | fundemental)
% b - vector
% xnorm - L2 norm minimization outcome (shold be a real positive scalar)
% x - solution
%
% example:

% Data
n=2;
m=20;
A=rand(m, n);
b=rand(m, 1);
xnorm = 2;

% solution
xopt = quadMin(A, b, xnorm);

% Check
theta=linspace(-pi, pi, 100);
% Generate 100 vectors such that |x|=xnorm
X=xnorm * [cos(theta); sin(theta)];
R=bsxfun(@minus, A * X, b);
r=sum(R.^2, 1); % r = J(x)

% Graphic
figure(1); clf
plot(theta, r);
hold on;
ttopt=atan2(xopt(2), xopt(1));
line([ttopt, ttopt], [0, max(r)], 'Color', 'r');
legend('cost function', 'optimal solution cost');
grid on;

% Dan I. Malta 2015
%}
function x = quadMin(A, b, xnorm)
    % housekeeping
    xnorm = real(xnorm);
    if xnorm == 0
        x = zeros(size(A,2), 1);
        return;
    end

    % formulate problem as an eigenvalue problem
    H = A' * A;
    g = A' * b;
    gn = g / xnorm;
    A2 = speye(size(H));
    A1 = -2 * H;
    A0 = H * H' - gn*gn';
    
    % eigenvalues
    [V, lambda] = polyeig(A0, A1, A2);  %#ok<ASGLU>
    
    % remove imaginary eigenvalue
    lambda = lambda(imag(lambda) == 0);
    if isempty(lambda)
        error('quadMin encountered problem with floating point accuracy'); %#ok<ERTAG>
    end

    % choose the smallest eigenvalue
    lambda = min(lambda);
    
    % solution
    x = (H - lambda * speye(size(H))) \ g;
end
