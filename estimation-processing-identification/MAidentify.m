%{
% MAidentify uses a direct and brute-force method to transfer function identification of moving average process.
% although specifically designed for MA process, works well for AR process (as can be seen in the example below)
%
% W = system output
% V = system input
% Ts = sampling perid
% n = estimated system order (number of poles if an all poles system assumed)    
% num, den = identified system transfer function numerator & denominator
%
% see ARidentify for adaptive identification of AR process.
%
% example:

% createa an AR(2) process whos input is a sine
Ts = 0.05;
sigma = 0.5;
t = 0 : Ts : 10;
y = sin(t);
for i = 3 : length(t)
    y(i) = 0.3 * y(i - 1) + 0.4 * y(i-2) + sigma * randn();
end

% estimate model
[num, den] = MAidentify(y, sin(t), Ts, 2);
out = filter(num, den, sin(t)) / num(1);

figure;
subplot(2,2,1);
plot(t, sin(t));
grid on;
title('model input');

subplot(2,2,2);
plot(t, y);
grid on;
title('model output (model IS AR(2) = 1-0.3*B-0.4*B^2 + N(0, \sigma^2=0.5))')

subplot(2,2,[3 4]);
y = sin(t);
for i = 3 : length(t)
    y(i) = y(i) + 0.3 * y(i - 1) + 0.4 * y(i-2);
end
plot(t, y, 'r', t, out, 'b');
grid on;
title('model without noise vs. estimated model output (estimating an AR process using a method for estimating MA process!!)');
legend('clean AR model (1-0.3*B-0.4*B^2) output', 'estimated model output');

%    
% Dan I. Malta 2011    
%}
function [num, den] = MAidentify(W, V, Ts, n)
    % recording duration
    T = length(W) * Ts;

    % time window length
    Tw = (T - (2 * n + 1) * Ts);
    
    % number of time window samples
    m = round(Tw / Ts);

    % estimate a unified process & sensor noise model
    for i = m : -1 : 1
        for j = 1 : n
            Q(m - i + 1, j) = W((n+1) + i - j);
        end
        Vk(m - i + 1) = V(i + (n + 1));
        Wk(m - i + 1) = W(i + (n + 1));
    end

    % estimated state space system
    X = [Q Vk'];
    A = (X' * X)^-1 * X' * Wk';
    A1 = A';
    
    % estimated numerator & denominator
    num = [A1(end) zeros(1:n)];
    den = [1 -A1(1:n)];
end
