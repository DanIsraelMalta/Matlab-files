%{
% fftSmooth performs a fast and efficient robust smoothing by minimizing
% the general cross validation function, in the frequency response!
%
% x - vector to be smoothed
% y - smoothed x
%
% example:

t = linspace(0, 2 * pi, 1000);
x = 2 * cos(t) .* (1-cos(t)) + 0.7 * randn(size(t)) ;
y = fftSmooth(x);
figure;
plot(t, x, 'r.', t, y,'k','linewidth',2);
grid on;

%
% Dan I. Malta 2010
%}
function y = fftSmooth(x)
    % housekeeping
    n = length(x);
    siz0 = size(x);
    x = x(:).';
    
    % transform to the frequency domain
    Lambda = 2 - 2 * cos(2 * (0 : n - 1) * pi / n);
    Y = fft(x);
    
    % smoothing optimization
    fminbnd(@GCVscore, -10, 30, optimset('TOLX',.1));
    
    % return tot ime domain
    if isreal(Y)
        y = ifft(Gamma .* Y, 'symmetric');
    else
        y = ifft(Gamma.*Y);
    end
    
    % output reshape
    y = reshape(y, siz0);

    % general cross validation function
    function GCVs = GCVscore(xi_p)
        s = 10 ^ xi_p;
        Gamma = 1 ./ (1 + s * Lambda.^2);
        RSS = norm(Y .* (Gamma - 1))^2;
        TrH = sum(Gamma);
        GCVs = RSS / (1 - TrH / n)^2;
    end
end
