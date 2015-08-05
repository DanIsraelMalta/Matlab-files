%{
% SDF is a variable mass rigid body six degrees of freedom equation of motion solver.
% SDF uses ODE87 for improved acuracy (thus running time is high).
%
% "mandatory" inputs:
% - xi_focrces - 3x1 vector of forces (BODY)
% - xi_moments - 3x1 vectory of moments (BODY)
% - xi_dmass - mass rate of change for the aircraft
% - xi_dinertia - 3x3 inertia tensor matrix rate of change
% - xi_time - time series vector
% "optional" inputs:
% - initial position - 3x1 vector of initial position (assumed [0, 0, 0])
% - initial velocity - 3x1 vector of initial velocity (BODY) (assumed [0, 0, 0])
% - initial body rates - 3x1 vector of initial rates (BODY) (assumed [0, 0, 0])
% - initial mass - aircraft initial mass (assumed 430kg)
% - initial inertia - 3x3 inetia matrix
% - initial attitude - 3x1 vector of initial euler angles (assumed [0, 0, 0])
% output:
% - xo_t - simulation time
% - xo_y - 1x19 output vector, as follows:
%                 y(1:3) = rates
%                 y(4:6) = vel
%                 y(7:9) = pos
%                 y(10) = Mass
%                 y(11:19) = Inertia
%
% usage example:

% rigid body
 XYZ = [0 0 0]; % initial position
 UVW = [50 0 -5]; % initial velocities (body)
 PsiThetaPhi = [0 45 0] * pi / 180; % initial attitude (euler)
 PQR = [0 5 0] * pi / 180; % initial rates (body)
 M = 10.0; % inital mass
 I = eye(3); % initial inertia
 LMN = [0 0 0]'; % acting moments
 F_XYZ = [-2 0 -9]'; % acting frces (body)
 dM = -0.01; % mass rate of loss
 dI = -0.01 * ones(3); % inertia rate of loss

% simulation
 tout = [0, 6.3];
 [t_out, y_out]=SDF(F_XYZ, LMN, dM, dI, tout, XYZ, UVW, PQR, M, I, PsiThetaPhi);

% visualization
figure;
subplot(2,2,1);
plot(t_out, y_out(:, 4));
ylabel('U [m/s]');
title('A not-so-perfect boomerang throw');
grid on;
subplot(2,2,2);
plot(t_out, y_out(:, 6));
ylabel('W [m/s]]');
grid on;
subplot(2,2,3);
plot(t_out, y_out(:,9));
ylabel('Ze [m]');
grid on;
subplot(2,2,4);
plot(t_out, y_out(:,7));
ylabel('Xe [m]');
grid on;

%
% Dan I. Malta 2011
%}
function [t,y] = SDF(xi_forces, xi_moments, xi_dmass, xi_dinertia, xi_time, varargin )
if nargin < 5
    error('Less than the minimum 5 inputs.') %#ok<ERTAG>
elseif nargin > 11
    error('More than maximum 11 inputs.') %#ok<ERTAG>
else
    % initialize optional inputs
    initialPosition = [0 0 0];
    initialVelocityB = [0 0 0];
    initialRatesB = [0 0 0];
    initialMass = 430.0;
    initialInertia = eye(3);
    initialEuler = [0 0 0]; %#ok<NASGU>
    
    % check optional inputs
    switch nargin
        case 6
            initialPosition = varargin{1};
        case 7
            initialPosition = varargin{1};
            initialVelocityB = varargin{2};
        case 8
            initialPosition = varargin{1};
            initialVelocityB = varargin{2};
            initialRatesB = varargin{3};
        case 9
            initialPosition = varargin{1};
            initialVelocityB = varargin{2};
            initialRatesB = varargin{3};
            initialMass = varargin{4};
        case 10
            initialPosition = varargin{1};
            initialVelocityB = varargin{2};
            initialRatesB = varargin{3};
            initialMass = varargin{4};
            initialInertia = varargin{5};
        case 11
            initialPosition = varargin{1};
            initialVelocityB = varargin{2};
            initialRatesB = varargin{3};
            initialMass = varargin{4};
            initialInertia = varargin{5};
            initialEuler = varargin{6};
    end
    
    % DCM construction
    cosPhi = cos(initialEuler(1));       sinPhi = sin(initialEuler(1));
    cosTheta = cos(initialEuler(2));  sinTheta = sin(initialEuler(2));
    cosPsi = cos(initialEuler(3));       sinPsi = sin(initialEuler(3));
    dcm = [cosTheta*cosPsi,                                                              cosTheta*sinPsi,                                                        -sinTheta;...
                  sinPhi*sinTheta*cosPsi-cosPhi*sinPsi,         sinPhi*sinTheta*sinPsi+cosPhi*cosPsi,  sinPhi*cosTheta; ...
                  cosPhi*sinTheta*cosPsi+sinPhi*sinPsi,        cosPhi*sinTheta*sinPsi-sinPhi*cosPsi,   cosPhi*cosTheta];

    % convert initialPosition into body coordinates
    initialPosition = (dcm * initialPosition')';
end

if initialMass <= 0
    error('Mass should be a positive non-zero value (this is not nuclear physics...).')     %#ok<ERTAG>
else
    y0 = [initialRatesB, initialVelocityB, initialPosition, initialMass, reshape(initialInertia, 1, 9)];
    options = odeset('MaxStep', 1/50, 'Mass', @massMat, 'MStateDependence', 'strong');
    [t,y] = ODE87(@DyDt, xi_time, y0, options, xi_forces, xi_dmass, xi_moments, xi_dinertia);
end

% calculate equations of motion rate
    function xo_dydt = DyDt(t, y, xi_forces, xi_dmass, xi_moments, xi_dinertia) %#ok<INUSL>
        % xo_dydt(1:3) = rates_dot
        % xo_dydt(4:6) = accel
        % xo_dydt(7:9) = vel
        % xo_dydt(10) = xi_dmass
        % xo_dydt(11:19) = xi_dinertia
        % y(1:3) = rates
        % y(4:6) = vel
        % y(7:9) = pos
        % y(10) = Mass
        % y(11:19) = Inertia

        xo_dydt = [xi_moments - cross(y(1:3), reshape(y(11:19),3,3) * y(1:3)) - xi_dinertia*y(1:3); ...
                                xi_forces - y(10)*cross(y(1:3),y(4:6)) - xi_dmass*y(4:6); ...
                               y(4:6);...
                               xi_dmass;...
                               reshape(xi_dinertia,9,1)];
    end

% calculate mass matrix
    function M = massMat(t, y, xi_forces, xi_dmass, xi_moments, xi_dinertia) %#ok<INUSD,INUSL>
        % initialized mass matrix
        M = ones(19,19);

        % fill matrix with state
        M(1,1) = y(11);
        M(1,2) = y(14);
        M(1,3) = y(17);
        M(2,1) = y(12);
        M(2,2) = y(15);
        M(2,3) = y(18);
        M(3,1) = y(13);
        M(3,2) = y(16);
        M(3,3) = y(19);
        M(4,4) = y(10);
        M(5,5) = y(10);
        M(6,6) = y(10);
    end
end
