%{
% A library of usefull "one-liners" accumelated over the years.
% several of these "one-liners" use my other matlab functions.
%
% Dan I. Malta
%}

% condense a given vector x by a factor n.
y_condense = @(x, n) (mean(reshape(x(1 : n * floor(length(x) / n)), n, floor(length(x) / n))));

% return the index of an element at vector x that is closest to  value val.
y_index = @(x, val) (find((abs(x - val) - min(abs(x - val))) == 0));
  
% first order centeral divided difference (aka "slope").
y_diff1 = @(x, y) (pchip(0.5 * (x(1 : end - 1) + x(2 : end)), diff(y) ./ diff(x), x));

% second order centeral divided difference (aka "curvature").
% uses y_diff1.
y_diff2 = @(x, y) (pchip(0.5 * (x(1 : end - 1) + x(2 : end)), diff(y_diff1(x, y)) ./ diff(x), x));

% filter a vector of angles (given in radians [-pi, pi]) such that wrap-around values are taken care of.
%  use fftSmooth.m for smoothing purposes.
y_smooth_angle = @(x) (atan2(fftSmooth(sin(x)), fftSmooth(cos(x))));

% zero-phase filter function (api equal to filter function)
y_filtftil = @(b, a, x) (fliplr(filter(b, a, fliplr(filter(b, a, z)))));

% return a given vector x after implementing upon it a "time domain hanning filter" of length n [samples].
y_hann = @(x, n) (y_filtftil((1 - cos(2 * pi * (1 : n) / n)) / sum(1 - cos(2 * pi * (1 : n) / n)), 1, x));

% exponential moving average (autorecrusive, 1st order FIR).
%  x = observations, ts = sampling period, tau = filter time constant
y_ema = @(x, ts, tau) (filter(ts / (ts + tau), [1, -tau / (ts + tau)], x));
y_ema_zp = @(x, ts, tau) (y_filtftil(ts / (ts + tau), [1, -tau / (ts + tau)], x));

% one zero one pole IIR low pass filter.
%  x = observations, ts = sampling period, tau = filter time constant
y_lpf = @(x, ts, tau) (filter([ts, ts], [ts + 2 * tau, ts - 2 * tau], x));
y_lpf_zp = @(x, ts, tau) (y_filtftil([ts, ts], [ts + 2 * tau, ts - 2 * tau], x));

% one zero one pole IIR high pass filter.
%  x = observations, ts = sampling period, tau = filter time constant
y_hpf = @(x, ts, tau) (filter([2 * tau ^ 2, -2 * tau ^ 2], [ts + 2 * tau, ts - 2 * tau], x));
y_hpf_zp = @(x, ts, tau) (y_filtftil([2 * tau ^ 2, -2 * tau ^ 2], [ts + 2 * tau, ts - 2 * tau], x));
