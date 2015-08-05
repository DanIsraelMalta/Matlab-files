%{
% ARidentify estiamtes an autoregresive process using a recrusive maximum
% likelihood estimation schema
%
% x - time series data
% p - estimated model order
% a - vector holding AR process parameters
%
% see MAidentify for adaptive identification of MA process.
%
% example:

% createa an AR(2) process whos input is a sine
Ts = 0.05;
sigma = 0.5;
t = 0 : Ts : 10;
y = sin(t);
for i = 3 : length(t)
    y(i) = 0.3 * y(i - 1) + 0.4 * y(i-2) + sigma * randn();
end

% estimate model
a = ARidentify(y, 2);
out = filter(1, a, sin(t));

figure;
subplot(2,2,1);
plot(t, sin(t));
grid on;
title('model input');

subplot(2,2,2);
plot(t, y);
grid on;
title('model output (model = 1-0.3*B-0.4*B^2 + N(0, \sigma^2=0.5))')

subplot(2,2,[3 4]);
y = sin(t);
for i = 3 : length(t)
    y(i) = y(i) + 0.3 * y(i - 1) + 0.4 * y(i-2);
end
plot(t, y, 'r', t, out, 'b');
grid on;
title('model without noise vs. estimated model output');
legend('clean model (1-0.3*B-0.4*B^2) output', 'estimated model output');

%
% Dan I. Malta 2011
%}
function a = ARidentify(x, p)
    % scaling
    x = x * 1e-6;
    
    % housekeeping
    N = length(x);
    S = zeros(p + 1, p + 1);
    a_aux = zeros(p + 1,p);
    a_aux(1, :) = 1;
    b_aux = ones(p + 1,p);
    e_aux = zeros(p, 1);
    p_aux = zeros(p, 1);
    MLE = zeros(3, 1);
    pos = 1;

    % S creation
    for i = 0 : p
        for j = 0 : p      
            for n = 0 : N - 1 - i - j
                S(i + 1, j + 1) = S(i + 1, j + 1) + x(n + 1 + i) * x(n + 1 + j);
            end
        end
    end

    % a_uax & b_aux calculation
    e0 = S(1, 1);
    c1 = S(1, 2);
    d1 = S(2, 2);
    coef3 = 1;
    coef2 = ((N - 2) * c1) / ((N - 1) * d1);
    coef1 = -(e0 + N * d1) / ((N - 1) *d1);
    ti = -(N * c1) / ((N - 1) * d1);
    raices = roots([coef3 coef2 coef1 ti]);
    for o=1:3
        if raices(o) > -1 && raices(o) < 1
            a_aux(2, 1) = raices(o);
            b_aux(p + 1, 1) = raices(o);
        end
    end
    e_aux(1, 1) = S(1,1)+2 * a_aux(2, 1) *S(1, 2) + (a_aux(2, 1)^2) * S(2, 2);
    p_aux(1, 1) = e_aux(1, 1) / N;

    % coefficient calculation
    for k=2:p
        Ck = S(1 : k, 2 : k + 1);
        Dk = S(2 : k + 1, 2 : k + 1);
        ck = a_aux(1:k,k-1)' * Ck * b_aux(p+1 : -1 : p+2-k, k-1);
        dk = b_aux(p+1 : -1 : p+2-k, k-1)' * Dk * b_aux(p+1 : -1 : p+2-k, k-1);
        coef3re = 1;
        coef2re = ((N - 2*k) * ck) / ((N-k) * dk);
        coef1re = -(k * e_aux(k-1, 1) + N * dk) / ((N-k) * dk);
        tire = -(N * ck) / ((N-k) * dk);
        raices = roots([coef3re coef2re coef1re tire]);
        for o = 1:3
            if raices(o,1) > -1 && raices(o,1) < 1
                MLE(o, 1) = ((1 - raices(o)^2) ^ (k/2)) / (((e_aux(k-1) + 2*ck*raices(o) + dk * (raices(o)^2))/N) ^ (N/2));
            end
        end
        [C, I] = max(MLE); %#ok<ASGLU>
        k_max = raices(I);
        for i = 1 : k - 1
            a_aux(i + 1, k) = a_aux(i + 1, k - 1) + k_max * a_aux(k - i + 1, k - 1);
        end
        a_aux(k + 1, k) = k_max;
        b_aux(p + 1 - k : p + 1, k) = a_aux(1 : k + 1, k);
        e_aux(k, 1) = e_aux(k - 1, 1) + 2 * ck * k_max + dk * k_max^2;
        p_aux(k, 1) = e_aux(k, 1) / N;
    end
    
    a = a_aux(:, p)';
end
