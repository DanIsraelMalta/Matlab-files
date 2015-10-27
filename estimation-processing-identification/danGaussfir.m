%{
% danGaussfir returns the coeeficients of gaussian filter.
% (same as Signal Processing Toolbox - gaussfir command)
%
% bandwidth - 3db bandwidth [rad/sec]
% oversampling - number of samples per period
% n - periods between start of filter impulse response and its peak
% Ts - sampling period
% h - filter coefficens array (normalized such that their sum is 1), whos length is 2 * oversampling * n + 1
%
% example:

% noisy signal
Ts = 0.01;
x = 0 : Ts : 4;
y = sin(pi * x) + cos(3 * pi * x) + sin(5 * pi * x.^2) + 0.2*randn(size(x));

% filter
h = danGaussfir(0.5, 2, 4, Ts);
ys = filter(h, 1, y);

% output
figure;

subplot(2,1,1);
plot(x, y, 'r', x, ys, 'k');
grid on;
legend('observation', 'filtered');
xlabel('Time');

subplot(2,1,2);
L = numel(x);
NFFT = 2^nextpow2(L);
Y = fft(y, NFFT) / L;
Ys = fft(ys, NFFT) / L;
f = linspace(0, 1, NFFT / 2) / 2 / Ts;
plot(f,2*abs(Y(1:NFFT/2)), 'r', f,2*abs(Ys(1:NFFT/2)), 'k');
grid on;
legend('observation', 'filtered');
xlabel('Frequency');

% Dan I. Malta 2015
%}
function h = danGaussfir(bandwidth, oversampling, n, Ts)
    delta = sqrt(log(2)) / (2 * pi * bandwidth * Ts);
    t = linspace(-n, n, oversampling * n);
    h = exp(-t.^2 / (2 * delta^2)) / (delta * sqrt(2 * pi));
    h = h / sum(h);
end
