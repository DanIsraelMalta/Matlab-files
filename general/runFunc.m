%{
% runFunc employs a general function (1 or 2 parameters) along a given
% window size on a given vector.
%
% xi - vector
% win - window length (number of points) for func
% func - inline (string) function to be implemented
%               should operate on a variable 'x'
% direction - window for running the function:
%                    'back' - run the function from current point backward
%                    'center' - run the function such that current point is
%                    in the center of the window
% backward - if '1', then run the function forward and then backward
% y - output
%
% example:

t = 0 : 0.01 : 4 * pi;
y = sin(t) + 0.1 * randn(size(t));
y_rms = runFunc(y, 20, 'norm(x) / sqrt(length(x))', 'center', 1);
y_avg3 = runFunc(y, 20, 'median(x)', 'center', 1);
y_max = runFunc(y, 10, 'max(x)',  'back', 1);
y_min = runFunc(y, 10, 'min(x)',  'back', 1);


figure;
plot(t, y, 'r', t, y_rms, 'b', t, y_avg3, 'k', t, y_max, 'g', t, y_min, 'g');
legend('signal', 'signal RMS (20 taps)', 'signal median filtered (20 taps)', 'signal min/max envelope');
grid on;

%
% Dan I. Malta 2008
%}
function yo = runFunc(xi, win, func, direction, backward)
    % housekeeping
    n = length(xi);
    f = inline(func);
     yo = xi;
    
    % backward window
    if strcmp(direction, 'back')
        % forward iteration
        for i = win + 1 : n
            yo(i) = f(xi(i : -1 : i - win));
        end

        % backward iteration
        if backward
            for i = n - win - 1 : -1 : 1
                yo(i) = f(xi(i : i + win));
            end
        end
    elseif strcmp(direction, 'center')  %center window
        % half window
        wlen = win / 2;
        if rem(win, 2) > 0 % odd length window
            wlen = wlen - 1;
        end
        
        % forward iteration
        for i = 1 : n
            yo(i) = f(xi(max(1, i - wlen) : 1 : min(n, i + wlen)));
        end

        % backward iteration
        if backward
            for i = n : -1 : 1
                yo(i) = f(xi(max(1, i - wlen) : 1 : min(n, i + wlen)));
            end
        end
    end
end
