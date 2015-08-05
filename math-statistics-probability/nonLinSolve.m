%{
%nonLinSolve uses the Jacobian-Free Newton-Krylov method for a set of non-linear equations
%
% F - function handle of nonlinear equations, in their residual form
% x0 - initial guess
% epsilon - tolerance
% max_iter - maximum number of iterations
% x - solution
% R - residuals
%
% example:

 F = @(x) [...
           % equation 1
           0.1309 * 1 / sqrt(3) * (1 / sqrt(1 / 3 + x(1) ^ 2 + x(1) ^ 2) + 1 / sqrt(1 / 3 + x(2) ^ 2 + x(4) ^ 2) +  1 / sqrt(1 / 3 + x(3) ^ 2 + x(3) ^ 2) + 1 / sqrt(1 / 3 + x(4) ^ 2 + x(2) ^ 2)) - 0.4352;
           % equation 2
           0.1309 * (x(1) / sqrt(1 / 3 + x(1) ^ 2 + x(1) ^ 2) + x(2) / sqrt(1 / 3 + x(2) ^ 2 + x(4) ^ 2) + x(3) / sqrt(1 / 3 + x(3) ^ 2 + x(3) ^ 2) + x(4) / sqrt(1 / 3 + x(4) ^ 2 + x(2) ^ 2)) - 0.1751;
           % equation 3
           0.1309 * (x(1) / sqrt(1 / 3 + x(1) ^ 2 + x(1) ^ 2) + x(4) / sqrt(1 / 3 + x(2) ^ 2 + x(4) ^ 2) + x(3) / sqrt(1 / 3 + x(3) ^ 2 + x(3) ^ 2) + x(2) / sqrt(1 / 3 + x(4) ^ 2 + x(2) ^ 2)) - 0.1751;
           % equation 4
           0.1309 * 1 / sqrt(3) * (x(1) / (1 / 3 + x(1) ^ 2 + x(1) ^ 2) + x(2) / (1 / 3 + x(2) ^ 2 + x(4) ^ 2) + x(3) / (1 / 3 + x(3) ^ 2 + x(3) ^ 2) + x(4) / (1 / 3 + x(4) ^ 2 + x(2) ^ 2)) - 0.1395;
          ];
 x0 = [0.4330, 0.4330, 0.1443, 0.1443];
 epsilon = 1e-5;
 max_iter = 10;
[x, R] = nonLinSolve(F, x0, epsilon, max_iter);

%
% Dan I. Malta 2014
%}
function [x, R] = nonLinSolve(F, x0, epsilon, max_iter)
    R = F(x0);         % initial residual
    x = x0;              % initialize solution vector
    counter = 1; % iteration counter
    
    while norm(R, 2) > epsilon
        j_v_approx = @(v)jacobApprox(v, F, x);
        v = gmres(j_v_approx, R);                             % solve for Krylov vector
        x = x - v';                                                                 % updated solution
        R = F(x);                                                                    % new residual
        
        if counter > max_iter
            error('nonLinSolve method did not converge! try increasing number of iterations...'); %#ok<ERTAG>
        end
        
        counter = counter + 1; % update iteration counter
    end

    fprintf('Error norm is: %f', norm(R, 2)); % norm of the residual
end

% "jacobian" pertubation
function y = jacobApprox(v, F, x)
    dim = size(x, 2);
    
    % calculate perturbation
    if norm(v, 2) > eps
        sum = 0;

        for i = 1 : dim
            sum = sum + sqrt(eps) * (1 + x(i));
        end

        per = (1 / (dim * norm(v, 2))) * sum;
    else
        sum = 0;

        for i = 1 : dim
            sum = sum + sqrt(eps) * (1 + x(i));
        end

        per = sum / dim;
    end

    R = F(x);                                % unperturbed residual
    xper = x' + per * v;         % perturbed vector
    Rper = F(xper);                % perturbed residual
    y = (Rper - R) / per;          % approximation of jacobian action on krylov vector
end
