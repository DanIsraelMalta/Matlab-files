%{
% RTSsmoothe performs the rauch-tung-striebel smoothing scehma.
% This is actualy a two stage (forward & backward) kalman smoother..
%
% x - signal to be smoothed
% Q - process noise (sigma square)
% R - measurement noise (sigma square)
% N - number of samples for second stage backward smoothing
%       if not specified, then taken as the length of x
% y - smoothed x
%
% example:

% noises
R = 0.5;
Q = R / 10;

% model
t = 0 : 0.1 : 12*pi;
y = sin(t) ;
i = 2 : 1 : length(t);
y(i) = y(i-1) + sqrt(R) * randn(size(i));

% r-t-s smooth
z = RTSsmooth(y, Q, R);

% visualize
figure;
plot(t, y, 'r', t, sin(t), 'b', t, zt, 'k');
legend('noisy observation', 'signal', 'smoothed observation');
grid on;

%
% Dan I. Malta 2015
%}
function y = RTSsmooth(x, Q, R)
    % columnize
    x = x(:);
    
    % housekeeping
    N = length(x);
    Ppred = zeros(size(x));
    Pcor = zeros(size(x));
    ypred = zeros(size(x));
    y = zeros(size(x));
    
    % forward pass initial step (kalman filtering)
    Ppred(1) = var(x);
    K = Ppred(1) / (Ppred(1) + R);
    ypred(1) = x(1);
    y(1) = ypred(1) + K * (x(1) - ypred(1));
    Pcor(1) = (1 - K) * Ppred(1);
    
    % forward pass iteration
    for i = 2 : N
        Ppred(i) = Pcor(i - 1) + Q;
        ypred(i) = y(i - 1);
        K = Ppred(i) / (Ppred(i) + R);
        y(i) = ypred(i) + K * (x(i) - ypred(i));
        Pcor(i) = (1 - K) * Ppred(i);
    end
    
    % backward pass (smoothing)
    for i = N - 1 : -1 : 1
        A = Pcor(i) / Ppred(i + 1);
        y(i) = y(i) + A * (y(i + 1) - ypred(i + 1));
    end
end
