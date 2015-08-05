%{
% Impulse smoother (direct approach based upon first moment; no need for t/chi2 distribution)
%
% remark:
% this smoother is based upon a heuristic concept I developed and tried to incorporate to a RT code (rejected by safety committee).
% Its main disadvatage is signal discretization, therfor it is adviced to pass its output thru a lowpass filter.
%
% x      = observation
% win = window length [samples] for impulse smoothing (should be chosen as impulse window width + 2)
% y       = x after it was impulse smoothed
%
% example:

% noisy filter
ts = 0.001;
t = 0 : ts : 2*pi;
len = numel(t);
y = sin(t);
y(len / 10 : 30 : len - 10) = 1.0;
y(len / 15 : 45 : len - 15) = -1.0;

%impulse removal
win = 3;
z = impulseSmooth(y, win);

% discretization removal utilizing double-sided triangular lowpass filter
z_lpf = filter([1 : win, win - 1 : -1 : 1] / win^2, 1, z);
zl = fliplr(z_lpf);
z_lpf = filter([1 : win, win - 1 : -1 : 1] / win^2, 1, zl);
z_lpf = fliplr(z_lpf);

figure;
plot(t, y, 'r', t, z, 'k', t, z_lpf, 'b');
grid on;
legend('observation', 'impulse smooth', 'impulse smooth and the low pass filtered');

%
% Dan I. Malta 2014
%}
function y = impulseSmooth(x, win)
    % housekeeping
    y = x;
    len = numel(x);
    win2 = win^2;
    
    % impulse removal
    if (len > win)
        for i = win + 1 : win :len
            s = x(i - win : i);
            e = median(s);
            p = sum(s > e);
            n = sum(s < e);
            t = abs(sum(s(s > e)) - e);
            smooth = e + (p - n) * t / win2;
            y(i - win : i) = smooth;
        end
    end
end
