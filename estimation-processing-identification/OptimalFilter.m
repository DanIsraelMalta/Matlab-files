%{
% OptimalFilter returns a direct form optimal filter (of any order)
% coefficients for low pass FIR filters, based upon Durrani-Chapman method
%
% note that the original algorithm suffers from serious numerical problems, so I
% had to scale things up/down and change various calculations to accommodate them.
%
% order - filter order (specify 1 for IIR filter)
% Fc - filter cutoff point / bandwidth [Hz]
% Fs - sampling frequency [Hz]
%
% aa - filter numerator
% bb - filter denominator
% the function also plots the filter characteristics.
%
% example:

% design an optimal 2nd order filter with cutoff frequency of 0.8Hz
% and apply to on a noisy step input
clc;
Ts = 0.05;           % sampling period [sec]
fc = 0.8;             % cutoff frequency [hz] (frequencies beneth this value should be filtered out)
order = 2;        % filter order in poles
Step = [zeros(1,50), 0.5*ones(1,5), ones(1,100)];
Step = Step+0.05*randn(size(Step));
Time = (1:length(Step))*Ts;
t0 = cputime;
[aa bb] = OptimalFilter(order,fc,1/Ts);
t=cputime-t0;
disp(['Optima Filter was found in ', num2str(t), ' [sec].']);
disp(['Optimal Filter requirements: F_c = ', num2str(fc) ,', F_s = ', num2str(1/Ts), ', order = ', num2str(order), '.']);
disp(['Optimal Filter:Y[n] = (',num2str(1-bb),')*Y[n-1] + ',num2str(bb),'*X[n]']);
disp(['Optimal Filter Time Constant: tau = ', num2str(((1-bb)/bb)*Ts), ' [sec]']);
fil = filter(bb,aa,Step);
figure;
plot(Time, Step, 'k', Time, fil, 'b');
grid on; legend('Step Cmd+ Blimp (during step occurrence) + White Noise',...
                                    ['LPF (Y[n] = (',num2str(1-bb),')*Y[n-1] + ',num2str(bb),'*X[n])']);
xlabel('Time'); ylabel('Step Cmd');

%
% Dan I. Malta 2009
%}
function [aa bb] = OptimalFilter(Order, Fc, Fs)
    % housekeeping
    fnot = Fc / Fs;
    N = Order;
    format long; % the algorithm is "numerical accurate"
    
    % create the Discrete Prolate Spheroidal Sequences
    sigma = zeros(N, N);
    for k = 1 : N - 1
        sigma(k, k) = (N + 1 - 2 * k)^2 / 4 * cos(2 * pi * fnot);
        sigma(k, k + 1) = 0.5 * k * (N - 1 - (k - 1));
        sigma(k + 1, k) = 0.5 * k * (N - k);
    end
    sigma(N, N) = (N + 1 - 2 * N)^2 / 4 * cos(2 * pi * fnot);
    [v lambda] = eig(sigma);
    lambda = nonzeros(lambda);
    [theta ind] = min(lambda); %#ok<ASGLU>
    v = v(:, ind);
    v = v ./ v(1);
    
    % filter gain and argument
    f = fnot;
    R = zeros(N, N);
    for k = 1 : N
        for l = 1 : N %#ok<FORPF>
            R(k, l) = exp(i * 2 * pi * f * (k - l));
        end
    end
    Arg = v' * R * v;
    K = sqrt(1 / Arg);
    
    % create polynom
    denom = K^2 * conv(v, v);
    denom = denom / sign(denom(1)) - [zeros(N - 1, 1); ones(1, 1); zeros(N - 1, 1)];
    denom = denom ./ denom(1);
    fact = roots(denom);
    stablepoles = nonzeros(fact .* (abs(fact) < 1));
    
    % scale factor (numerical issues)
    f = (0 : 1000) / 1000;
    
    % locate constant gain
    H = 1;
    for k = 1 :length(stablepoles) %#ok<FORPF>
        H = (exp(i * 2 * pi * f) - stablepoles(k)) .* H;
    end
    
    % transfer function
    H = 1 ./ H;
    G = 1 / abs(H(1));
    H = G * H';
    
    % extract only the stable poles
    aa = poly(stablepoles);
    bb = G;
    
    % Plot the poles and zeros
    figure;
    subplot(2,1,1);
    plot(real(fact), imag(fact), 'x');
    hold on;
    t = 0 : 0.05 : (2 * pi);
    plot(cos(t), sin(t));
    plot(0, 0, 'o');
    grid on;
    axis([-1.2 1.5 -1.2 1.2]);
    xlabel('Real Part');
    ylabel('Imaginary Part');
    
    % plot magnitude and phase
    subplot(2,1,2);
    MagnResp = 20 * log10(abs(H));
    phase = (180 / pi) * unwrap(angle(H)); %phase = phase-phase(1);
    [AX,H1,H2] = plotyy((0 : 1000) * Fs / 1000, MagnResp, (0 : 1000) / 1000 * Fs, phase, 'plot'); %#ok<NASGU,NASGU>
    title('Magnitude & Phase Response Characteristics');
    set(get(AX(1),'Ylabel'),'String','Magnitude response [dB]');
    set(get(AX(2),'Ylabel'),'String','Phase [\circ]');
    set(get(AX(1),'Xlabel'),'String',['Frequency (Hz), F_c = ' num2str(Fc) ', and F_s = ' num2str(Fs)'.']);
    set(get(AX(2),'Xlabel'),'String',['Frequency (Hz), F_c = ' num2str(Fc) ', and F_s = ' num2str(Fs)'.']);
    set(AX(1), 'Xlim', [0 Fs/2]); set(AX(2), 'Xlim', [0 Fs / 2]);
    set(AX(1), 'XScale', 'log'); set(AX(2), 'XScale', 'log');
    grid on;
end
