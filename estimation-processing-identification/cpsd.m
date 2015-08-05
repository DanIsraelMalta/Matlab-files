%{
% cpsd calculates the cross power spectral density, coherence & phase of
% two input signals.
% cpsd uses wpsd.m.
%
% x, y - input signal
% window - window weights (should be the same length as nfft) [samples]
% noverlap - number of x samples that are common to two adjacent segments
% nfft - fourier transform length (should be the same length as window) [samples]
% fs - x sampling frequency [Hz]
% Pxx - x periodegram
% Pyy - y periodegram
% Pxy - cross periodegram
% Phase - cross phase
% C - coherence
% F - frequency vector
%
% Dan I. Malta 2011
%}
function [Pxx, Pyy, Pxy, Phase, C, F] = cpsd(x, y, window, noverlap, nfft, fs)
    [Pxx, F] = wpsd(x, window, noverlap, nfft, fs);
    [Pyy, F] = wpsd(y, window, noverlap, nfft, fs);
    Pxy = Pxx ./ Pyy;
    Phase = atan2(imag(Pxy), real(Pxy));
    C = Pxy .* conj(Pxy) ./ (Pxx .* Pyy);
end
