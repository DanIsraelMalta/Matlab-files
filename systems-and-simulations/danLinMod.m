%{
% danLinMod - linmod for any sort of non-simulink nonlinear dynamic system
%
% for a non linear dynamic system : {dx/dt = f(x, u), y = h(x, u)}
% or a nonlinear discrete system: {x[k+1] = f(x[k], u[k]), y[k] = h(x[k], u[k])}
% danLinMod(f, x,u ,h) will linearize th system at the reference point
% specified by the state vector (x) and input vector (u) and return the
% appropriate matrix corresponding to the linearized state space model:
% {dx/dt = Ax + Bu, y = Cx + Du} or {x[k+1] = Ax[x] + Bu[x], y[x] = Cx[x] + Du[x]}
% and {dx = f(x, u), y = h(x, u)}.
%
% linearization algorithm uses a complex step differentiation (since the
% jacobian is calculated for a non linear function).
%
% notes:
% 1) 2008: tested it against simulink models with dlinmod and kantor ssdf and it
%     gives correct results up to fourth digits after the point.
% 2) used this function for acceleration control loop design on H450.
%     it gave me bad results when compared to baruch ssdf trim on ground.
%     but when compared to yaniv new ground.for / wheel.for modules - the results were accurate.
% 3) 2010: when linearizing a "non-linear" phugoid model for velocity estimaotr design,
%     the results seem to embody some sort of constant bias in the sensor transformation matrix.
%     solved it by changing the numerical algorithms from my nJacob.m to a complex step jacobian.
% 
% example:

% nonlinear system dynamics:
%  / dX[1]dt = X[2]^2 + U[2]^2;
%  \ dX[2]dt = U[1]^3 - 2*zeta*wn * X[2]^2 - wn*wn * X[1];
% nonlinear system output:
% / Y[1] = 0.5 * sqrt(X[1])
% \ Y[2] = 1.5 * pow(X[2], 0.75)

% nonlinear system definition
Zeta = 0.4;
Wn = 31.4;
f = @(x, u) [x(2)^2 + u(2)^2; u(1)^3 - 2*Zeta*Wn*x(2)^2 - Wn*Wn*x(1)];
h = @(x) [0.5 * x(1)^0.5; 1.5 * x(2)^0.75];

% states and input
x = [3; 3];
u = [-2; 2];

% linearization
[A, B, C, D, dx, y] = danLinMod(f, x, u, h);

% nonlinear system linearization aresult around {x, u}
linss = ss(A, B, C, D);
ltiview(linss);

%
% Dan I. Malta 2010
%}
function [A, B, C, D, dx, y]=danLinMod(f, x, u, h)
    % state equation as function of input vector
    if nargin < 3 || isempty(u)   % input-free state equation
        [A, dx] = complexJacobian(f, x);
        B = [];
        nu = 0;
    else
        [A, dx] = complexJacobian(@(x)f(x, u), x);
        B = complexJacobian(@(u)f(x, u), u);
        nu = numel(u);
    end
       
    % full state
    if nargin < 4
        C = eye(numel(x));
        D = zeros(numel(x), nu);
        y = x;
    else
        try     %output equation with input
            [C, y] = complexJacobian(@(x)h(x, u), x);
            D = complexJacobian(@(u)h(x, u), u);
        catch %#ok<CTCH> %output equations without input
            [C, y] = complexJacobian(h, x);              
            D = [];
        end
    end
end

% complex step jacobian matrix numerical calculation
function [A,z]=complexJacobian(fun, x)
    % housekeeping
    n = numel(x);                     % size of independent
    m = numel(z);                     % size of dependent
    A = zeros(m, n);                   % allocate memory for the Jacobian matrix
    
    % initialization
    z = fun(x);
    step = n * eps;                        % differentiation step size
    
    % loop over all independent variables
    for k = 1 : n
        x1 = x;
        
        % increment k'th independent variable
        x1(k) = x1(k) + step * i;
        
        % complex step differentiation
        A(:, k) = imag(fun(x1)) / step;
    end
end
