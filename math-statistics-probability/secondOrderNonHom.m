%{
% numeric solution of:
%
% d^2x(t)                  dx(t)
% ----------- = p(t) *  -------- + q(t) * x(t) + r(t)    , x(a) = alpha, x(b) = beta
%   dt^2                       dt
%
% xi_p, xi_q, xi_r, xi_a, xi_b, xi_alpha, xi_beta - as defined in the problem definition
% xi_steps - number of steps
%
% example:

p = @(t) (0.1*t-1);
q = @(t) sqrt(t);
r = @(t) (0.5*t.^0.2);
[t, x] = secondOrderNonHom(p, q, r, 0, 10, 3, 0, 1000);

figure;
plot(t, x, 'b', t, feval(p, t), '--r', t, feval(q, t), '--g', t, feval(r, t), '--m', 0, 3, 'ok', 10, 0, 'ok');
xlabel('Time [sec]');
legend('X(t)', '0.1*t-1', 'sqrt(t)', '0.5*t^{0.2}', 'boundary condition', 'Location', 'Best');
title('X''''(t) = (0.1*t-1)*X''(t) + sqrt(t)*X(t) + 0.5*t^{0.2}');
grid on;

%
% Dan I. Malta 2016
%}
function [xo_t, xo_x] = secondOrderNonHom(xi_p, xi_q, xi_r, xi_a, xi_b, xi_alpha, xi_beta, xi_steps)
    % allocation
    xo_t = zeros(1, xi_steps + 1); % solution time vector
    xo_x = zeros(1, xi_steps - 1); % solution vector
    Vb = zeros(1, xi_steps - 1); % right hand side vector
    Vd = zeros(1, xi_steps - 1); % diagonal vector
    Vt = zeros(1, xi_steps - 1); % time vector
    Va = zeros(1, xi_steps - 2); % sub diagonal vector
    Vc = zeros(1, xi_steps - 2); % super diagonal vector

    % constants
    h = (xi_b - xi_a) / xi_steps;
    hSqr = h^2;
    hHalf = h / 2;

    % triangular system vector build
    for j = 1 : xi_steps - 1
        Vt(j) = xi_a + h * j;
        Vb(j) = -hSqr * feval(xi_r, Vt(j));
        Vd(j) = 2 + hSqr * feval(xi_q, Vt(j));
    end
    for j = 1 : xi_steps - 2
        Va(j) = -1 - hHalf * feval(xi_p, Vt(j + 1));
        Vc(j) = -1 + hHalf * feval(xi_p, Vt(j));
    end

    % boundary value correction for right hand side vector
    Vb(1)   = Vb(1)   + (1 + hHalf * feval(xi_p, Vt(1))) * xi_alpha;
    Vb(xi_steps-1) = Vb(xi_steps-1) + (1 - hHalf * feval(xi_p, Vt(xi_steps - 1))) * xi_beta;

    % triangular linear system solution
    xi_steps = length(Vb);
    for k = 2 : xi_steps - 1
        mult = Va(k-1) / Vd(k-1);
        Vd(k) = Vd(k) - mult * Vc(k - 1);
        Vb(k) = Vb(k) - mult * Vb(k - 1);
    end
    xo_x(xi_steps - 1) = Vb(xi_steps - 1) / Vd(xi_steps - 1);
    for k = (xi_steps - 2) : -1 : 1
        xo_x(k) = (Vb(k) - Vc(k) * xo_x(k + 1)) / Vd(k);
    end

    % output
    xo_t = [xi_a, Vt, xi_b];
    xo_x = [xi_alpha, xo_x, xi_beta];
end
