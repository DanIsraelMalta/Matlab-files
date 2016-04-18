%{
% chirp Z transform
%
% xi_x    - input vector
% xi_fs   - sampling frequency [Hz]
% xi_zoom - [zoom start frequency, zoom end frequency] [Hz]
% xi_len  - number of points in zoom
% xo_czt - chirp Z transform
% xo_w - frequency vector
%
% Dan I. Malta 2016
%}
function [xo_czt, xo_w] = chirpZ(xi_x, xi_fs, xi_zoom, xi_len)
    % columnize input
    [xrow, xcol] = size(xi_x);
    if xrow == 1
        xi_x = xi_x(:);
        [xrow, xcol] = size(xi_x);
    end
    
    % 'chirp' final point
    B = xi_len + xrow - 1;

    % power-of-two FFT length
    nfft = 2 ^ nextpow2(B);

    % Z-transform 'chirp' contour path points ratio
    k = ((1 - xrow) : max(xi_len - 1, xrow - 1)).';
    k = (k .^ 2) / 2;
    W = exp(-2 * pi * j * (xi_zoom(2) - xi_zoom(1)) / (xi_len * xi_fs));
    W = W .^ k;
    
    % chirp Z transform starting point
    nn = (0 : xrow - 1)';
    A = exp(-2 * pi * j * xi_zoom(2) / xi_fs);
    A = A .^ (-nn);
    A = A .* W(xrow + nn);
    
    % convolution
    one = ones(1, xcol);
    y = xi_x .* A(:, one);
    yConv = fft(y, nfft);
    chirp = fft(1 ./ W(1 : B), nfft);
    yConv = yConv .* chirp(:, one);
    
    % 'chirp'
    xo_czt = ifft(yConv);
    xo_czt = xo_czt(xrow : B, :) .* W(xrow : B, ones(1, xcol));
    xo_czt = flipud(xo_czt);
    
    % frequency vector
    cztLen = length(xo_czt);
    xo_w = xi_zoom(1) + (0 : cztLen - 1)' * (xi_zoom(2) - xi_zoom(1)) / cztLen;
end
