%{
% danImpulse returns the responce of a given system to an impulse for a given time span
% it is the control toolbox equivalent of impulse command, but it accepts
% function handles instead of control toolbox LTI objects.
%
% system - a single variable function handle which defines a system in the Laplace 
%                    domain (where the single variables is the Laplace variable 's')
% time_span - time duration for the impulse response simulation
% out - output
%
% example:

% sampling period [sec]
Ts = 1/20;

% servo (damping ratio = 0.7, natural frequency = 31.4 [rad/sec])
servo = @(s) (31.4^2 / (s.^2 + 2*0.7*31.4*s+31.4^2));

% zero order hold
zoh = @(s) ((-Ts/6*s + 1) / (Ts/3*s + 1));

% elevator to pitch angle transfer function
elevatorToPitch = @(s) (-(3.24*s.^2 + 2.2319*s + 0.0296) ./ (s.^4 + 2.0912*s.^3 + 5.7781*s.^2 + 0.0473*s + 0.0443));

% system (open loop)
OL = @(s) (servo(s) .* elevatorToPitch(s) .* zoh(s));

% system impulse reponse
t = 0 : Ts : 300;
OLout = danImpulse(OL, t);

% plot
figure;
plot(t, OLout);
title('system response to impulse');

%  Dan I. Malta 2015
%}
function out = danImpulse(system, time_span)
    % housekeeping
    coeffLen = 20; % number of coefficient (must be even number)
    coeffLenHalf = coeffLen / 2;
    
    % calculate the Harald coefficients
    HaraldCoeff = zeros(1, length(coeffLen));
    for i = 1 : coeffLen
        z = 0.0;
        for num = floor((i + 1) / 2) : min(i, coeffLenHalf)
            z = z + ((num ^ (coeffLenHalf)) * factorial(2 * num)) / ...
                          (factorial(coeffLenHalf - num) * factorial(num) * factorial(num - 1) * ...
                          factorial(i - num) * factorial(2 * num - i));
        end
        HaraldCoeff(i) = (-1) ^ (i + coeffLenHalf) * z;
    end
    
    % inverse Laplace transform
    out = zeros(1, length(time_span));
    for i = 1 : length(time_span)
        sum = 0.0;
        ln2OverTime = log(2.0) / time_span(i);
        for num = 1 : coeffLen
            sum = sum + HaraldCoeff(num) * feval(system, ln2OverTime * num);
        end
        out(i) = sum * ln2OverTime;
    end
end
