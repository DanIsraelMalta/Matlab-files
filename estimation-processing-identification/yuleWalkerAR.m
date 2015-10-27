%{
% yuleWalkerAR calculate the auto-regressive model parameters of a given
% uniformaly sampled time series using theYule-Walker method.
%
% see ARidentify for another method of autoregressive model estimation.
%
% x - time series to be modeled
% n - autoregressive model order
% A - autoregressive mode parameters
%      whice means that output is of the form:
%      y(t) + A(1) * y(t-1) + A(2)*y(t-2) + ... + A(n) * y(t-n) = e(t), where e is white noise
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
a = yuleWalkerAR(y, 2)
out = filter(1, a, sin(t));

figure;
subplot(2,2,1);
plot(t, sin(t));
grid on;
title('model input');

subplot(2,2,2);
plot(t, y);
grid on;
title('model output (model = 1-0.3*B-0.4*B^2 + N(0, \sigma^2=0.5))')

subplot(2,2,[3 4]);
y = sin(t);
for i = 3 : length(t)
    y(i) = y(i) + 0.3 * y(i - 1) + 0.4 * y(i-2);
end
plot(t, y, 'r', t, out, 'b');
grid on;
title('model without noise vs. estimated model output');
legend('clean model (1-0.3*B-0.4*B^2) output', 'estimated model output');

% Dan I. Malta 2105
%}
function A = yuleWalkerAR(x, n)
    % columnize
    x = x(:);
    
    % n-lags autocorrelation matrix
    nFFT = 2 ^ (nextpow2(length(x)) + 1);  % FFT length
    F = fft(x-mean(x) , nFFT);                         % detrended FFT
    F = F .* conj(F);                                             % power spectrum
    R  =  ifft(F);                                                   % autocorrelation
    R  =  R(1:(n + 1));                                               % Retain non-negative lags
    R  =  R ./ R(1);                                                     % Normalize autocorrelation
    R  =  real(R);                                                    % real valued autocorrelation
    
    % left side of Yule-Wlaker equation (autocorrelation in toeplitz form)
    YWleft = toeplitz(R(1:n));
    
    % right side of Yule-Wlaker equation
    YWright= -R(2 : n + 1);
    
    % autoregressive parameters
    A = [1 (YWleft \ YWright).'];
end
