%{
% bandPass performs perfect band pass filtering on a given time domain series
%
% {t, x} - time series
% {fmin, fmax} - frequeny band to be passed [Hz]
% y - band passed signal
%
% example:

% 50hz & 70hz sinosuid
Fs = 1000;
T = 1 / Fs;
L = 1000;
t = (0:L-1) * T;
x = 0.7*sin(2*pi*50*t) + 0.9*sin(2*pi*70*t) + 0.3 * randn(size(t));

% retrieve only the 50hz sinosuid
y = bandPass(t, x, 45, 55);

% vsualization
figure;
subplot(2,1,1);
plot(t, x, 'r', t, y, 'b');
legend('original signal', 'band passed signal');
title('Time Domain');
grid on;
subplot(2,1,2);
f = (Fs / 2) * linspace(0, 1, 128);
yf = fft(y, 256);
xf = fft(x, 256);
plot(f, 2 * abs(xf(1 : 128)), 'r', f, 2 * abs(yf(1:128)), 'b');
legend('original signal', 'band passed signal');
ylabel('magnitude'); xlabel('frequency [Hz]');
title('Frequency Domain');
grid on;

%
% Dan I. Malta 2007
%}
function y = bandPass(t, x, fmin, fmax)
    % housekeeping
    t = t(:);
    x = x(:);
    y = repmat(nan, size(x));
    
    % interpolate NaN's
    i = find(all(~isnan(x), 2));
    x = interp1(t(i), x(i,:), t(:));
    i = find(all(~isnan(x), 2));
    
    % frequency vector
    f = (0 : length(i) / 2) / (t(i(end)) - t(i(1)));
    
    % locate region of interest
    j = find(fmin <= f & f <= fmax);
    k = length(i) - j + 2;
    k(k > length(i)) = [];
    
    % filtering
    filt = zeros(length(i), size(x, 2));
    filt([j, k], :) = 1;
    if fmin > 0
        x(i, :) = detrend(x(i, :));
    end
    y(i, :) = real(ifft(filt .* fft(x(i, :))));
end
