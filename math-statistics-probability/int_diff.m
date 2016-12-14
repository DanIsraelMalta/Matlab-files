%{
% int_diff differentiates or integrates a signal, but in the frequency domain.
%
% xi_y    - signal to derive or integrate (time domain)
% xi_dt   - sampling interval
% xi_type - "diff" for derivative, "int" for integral
% xo_y    - outcome
%
% Dan I. Malta (2016)
%}
function xo_y = int_diff(xi_y, xi_dt, xi_type)
    % housekeeping
    len = length(xi_y);
    T   = xi_dt * len;
    xo_y = zeros(size(xi_y));

    % FFT (shifted)
    z  = fftshift(fft(xi_y));
    df = 1/T;

    % frequency grid
    if ~mod(len, 2)
        f = df * (-len / 2 : len / 2 - 1);
    else
        f = df * (-(len - 1) / 2 : (len - 1) / 2);
    end
    w = 2 * pi * f;

    % operation
    switch xi_type
        case 'diff'
            % inverse FFT (derivative by multiplying by i*omega)
            xo_y = ifft(ifftshift(z .* (1i * w)));
        case 'int'
            % Integrate in frequency domain by dividing spectrum with i*omega
            zn = z .* (-1i ./ w);
            
            % inverse FFT (make sure spectrum is conjugate symmetric!)
            xo_y = ifft(ifftshift(zn), 'symmetric');
        otherwise
            disp('Sorry, operation type is not recognized by int_diff.m.');
    end
end
