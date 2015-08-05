%{
% wpsd estimates a given signal power spectral density using Welch's method.
%
% x - input signal
% window - window weights (should be the same length as nfft) [samples]
% noverlap - number of x samples that are common to two adjacent segments
% nfft - fourier transform length (should be the same length as window) [samples]
% fs - x sampling frequency [Hz]
% Pxx - welch peridogram
% F - frequency vector
%
% Dan I. Malta 2010
%}
function [Pxx, F] = wpsd(x, window, noverlap, nfft, fs)
    % columnize
    x = x(:);
    window = window(:);

    % houskeeping
    len = length(x);
    df = fs / nfft;
    P = zeros(nfft, 1);
    xw = zeros(nfft, 1);
    
    % starting indices
    jump = nfft - noverlap;
    index = [1, jump : jump : len - nfft];
    if index(end) ~= len - nfft
        index = [index, len - nfft];
    end
    partitions = length(index);
    
    % welch
    for i = 1 : partitions
        xw = window .* x(index(i) : index(i) + nfft - 1);
        fourier = fft(xw, nfft);
        P = P +  fourier .* conj(fourier);
    end
        
    % periodogram
    Pxx = 2 * (P([2 2:nfft/2+1])) / (fs * partitions * norm(window)^2);
    F = (0 : max(size(Pxx)) - 1) * df;
end
