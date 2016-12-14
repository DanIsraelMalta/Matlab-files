%{
% secondDiscrete calculates the damping and natural frequency of a discrete second order system,
% given desired change in magnitude and pahse.
%
% xi_desired_phase     - desired phase
% xi_desired_magnitude - desired magnitude
% xi_dt                - sampling period [sec]
% xi_freq              - frequency vector [rad/sec]
% xo_zeta              - damping
% xo_wn                - natural ferquency
%
% Dan I. Malta (2016)
%}
function [xo_zeta, xo_wn] = secondOrderDiscrete(xi_desired_phase, xi_desired_magnitude, xi_dt, xi_freq)
    % negative phase
    if xi_desired_phase > 0
        xi_desired_phase     = -xi_desired_phase;
        xi_desired_magnitude = 1 / xi_desired_magnitude;
    end
    p = xi_desired_phase * pi / 180;
    
    % polynom coefficients
    c1 = cos(xi_dt * xi_freq);
    s1 = sin(xi_dt * xi_freq);
    c2 = cos(p + xi_dt * xi_freq);
    s2 = sin(p + xi_dt * xi_freq);
    c3 = cos(p + 2 * xi_dt * xi_freq);
    s3 = sin(p + 2 * xi_dt * xi_freq);
    c4 = cos(p);
    s4 = sin(p);
    
    % transfer function construction
    n = -(c3 - c1 / xi_desired_magnitude) * (s2 - s1 / xi_desired_magnitude) + (c2 - c1 / xi_desired_magnitude) * (s3 - s1 / xi_desired_magnitude);
    d = (c4 - c1 / xi_desired_magnitude) * (s2 - s1 / xi_desired_magnitude) - (c2 - c1 / xi_desired_magnitude) * (s4 - s1 / xi_desired_magnitude);
    a = -log(n / d) / 2 / xi_freq;
    cbt = (s3 - s1 / xi_desired_magnitude + (n / d) * (s4 - s1 / xi_desired_magnitude)) / (2 * sqrt(n / d) * (s2 - s1 / xi_desired_magnitude));
    b = acos(cbt) / xi_freq;
    
    % damping & frequency
    xo_zeta = sqrt(a^2 / (a^2 + b^2));
    xo_wn = b / sqrt(1 - xo_zeta^2);
    
    % complex number / zero damping handling
    if (imag(xo_zeta) ~= 0) | (imag(xo_wn) ~= 0) | (xo_zeta == 0) %#ok<OR2>
        xo_zeta = [];
        xo_wn   = [];
    else
        a = xo_zeta * xo_wn;
        b = xo_wn * sqrt(1 - xo_zeta^2);
        rts = roots([1 - 2 * exp(-a * xi_freq) * cos(b * xi_freq) exp(-2 * a * xi_freq)]);
        if any(abs(rts) >= 1)
            disp('non stable transfer function!');
        end
    end
end
