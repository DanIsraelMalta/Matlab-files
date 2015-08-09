%{
% TFidentify estimates the transfer function (numerator, denominator) of a given frequency response.
% TFidentify fit's a real (not complex) transfer function.
% TFidentify constraints the solution such that the denominator is stable.
%
% TFidentify requires stabilizePoly.m.
%
% H - system frequency response (could be complex)
% W - frequency vector on which H is given
%     to be "inlined" with various matlab functions, W is
%     given in [rad/sample] and within the interval [0, pi]
% NB - numerator order
% NA - denominator order
% Wt - vector of weights for specific frequencies (Wt length is equal to H & W length)
% ITER - allows the numerical estimation process maximal number of iterations
%        (minimal value, as taken by the dunction, is 30)
% TOL - stopage criteria for the norm of the gradient
%       (minimal value, as taken by the dunction, is 0.01)
% b - estimated transfer function numerator
% a - estimated transfer function denominator
%
% remarks:
% 1) if Wt is [], then optimization tries to minimize the |B-A*H|^2 cost function.
% 2) if Wt exists but ITER and TOL are [], then cost functio is Wt*|B-A*H|^2      ("error")
% 3) if Wt exists but ITER and TOL are not [], then cost function is Wt*|B/A-H|^2 ("output-error")
%
% TFidentify minimizes the cost function in the least-square sense, using two different algorithms:
% 1) in case of reamrk 2 - the identification algorithm is Levi's "complex curve fitting", it is performed
%    on the error equation, with my own modifications for real and stable transfer functions.
% 2) in case of remark 3 - the identification algorithm utilizes a simple iterative, nonlinear and uncosntrained
%    optimization, on the output-error equation.
%
% example:

% create a linear chirp
ts = 0.01;
t = 0 : ts : 9;
t0 = 0;
f0 = 0;
t1 = 1;
f1 = 1;
p = 1;
phi = 0;
beta   = (f1 - f0) .* (t1 .^ -1);
x = cos(2 * pi * ( beta./(1+p).*(t.^(1+p)) + f0.*t + phi/360)) + 0.5 * randn(size(t));

% system (impulse reponse coefficients)
h = [0.7, 0.5, 0.3, 0.1];

% output
y = filter(h, 1, x);
y = y + 0.25 * randn(size(y));

% direct frequency function calaulation
fs = 1 / ts;
Txy = fft(y, 256) ./ fft(x, 256);
Txy = Txy(1 : 128);
w = (fs / 2) * linspace(0, 1, 128);
w = w / (w(end) / pi);

% identify system
[b, a] = TFidentify(Txy, w, 3, 0, [], [], []);
yest = filter(b, a, x);

% plot
figure;

ax(1) = subplot(2, 2, 1);
plot(t, x);
title('system input (a noisy chirp)');
grid on;

subplot(2,2,2);
plot(1:length(h), h, 'or', 1:length(h), b, 'xb');
legend('system numerator coefficient', 'estimated numerator coefficient');
grid on;

ax(2) = subplot(2,2, [3 4]);
plot(t, y, 'r', t, yest, '--b');
legend('system output', 'estimation output');
grid on;

linkaxes(ax, 'x');

%
% Dan I. Malta 2012
%}
function [b,a]=TFidentify(H, W, NB, NA, Wt, ITER, TOL)
  % housekeeping
  if isempty(Wt)
    Wt = ones(size(W));
  end
  if isempty(ITER)
    ITER = 30;
  end
  if isempty(TOL)
    TOL = 0.01;
  end
  [rw, cw] = size(W);
  if rw > cw
    W=W';
  end
  [rg, cg] = size(H);
  if cg > rg
    H=H.';
  end
  [rwf, cwf] = size(Wt);
  if cwf > rwf
    Wt=Wt';
  end
  
  % initialization
  nk = 0;           % start with no zeros (all poles transfer function -> autoregressive process)
  nb = NB + nk + 1; % start with one pole excess
  Wt = sqrt(Wt);    % to fit the cost function

  % create exp(-j*omega)
  nm = max(NA, nb + nk - 1);
  OM = exp(-i * (0 : nm)' * W);
  
  % initial least squares solution (direct calculation)
  Dva = (OM(2 : NA + 1, :).') .* (H * ones(1, NA));
  Dvb = -(OM(nk + 1 : nk + nb, :).');
  D=[Dva Dvb] .* (Wt * ones(1, NA + nb));
  R = real(D' * D);                 % real transfer function
  Vd = real(D' * (-H .* Wt));       % real transfer function
  th = R \ Vd;
  
  % transfer function
  a = [1 th(1 : NA).'];
  b = [zeros(1, nk) th(NA + 1 : NA + nb).'];
  
  % Now for the iterative minimization
  indb = 1 : length(b);
  indg = 1 : length(a);
  
  % stabilize the denominator
  a = stabilizePoly(a);
  
  % initial curve fit
  GC = ((b * OM(indb, :)) ./ (a * OM(indg, :))).';
  e = (GC - H) .* Wt;
  Vcap = e' * e;
  t = [a(2 : NA + 1) b(nk + 1 : nk + nb)].';

  % minimization
  gndir = 2 * TOL + 1;
  l = 0;
  st = 0;
  while all([norm(gndir) > TOL, l < ITER, st ~= 1])
    l = l + 1;
    
    % gradient
    D31 = (OM(2 : NA + 1, :).') .* (-GC ./ ((a * OM(1 : NA + 1, :)).') * ones(1, NA));
    D32 = (OM(nk + 1 :nk + nb, :).') ./ ((a * OM(1 : NA + 1, :)).' * ones(1, nb));
    D3 = [D31 D32] .* (Wt * ones(1, NA + nb));
    
    % search direction (gauss-newton)
    e = (GC - H) .* Wt;
    R = real(D3' * D3);
    Vd = real(D3' * e);
    gndir = R \ Vd;
    
    % search along calculated direction (gndir; 20 steps)
    ll = 0;
    k = 1;
    V1 = Vcap + 1;
    while all([V1 > Vcap, ll < 20])
      t1 = t - k * gndir;
      if ll == 19
        t1=t;
      end
      
      % stable denominator
      a = stabilizePoly([1 t1(1 : NA).']);
      
      % numerator
      t1(1 : NA) = a(2 : NA + 1).';
      b = [zeros(1, nk) t1(NA + 1 : NA + nb).'];
      
      % dleast square and direction
      GC = ((b * OM(indb, :)) ./ (a * OM(indg, :))).';
      V1 = ((GC - H) .* Wt)' * ((GC - H) .* Wt);
      t1 = [a(2 : NA + 1) b(nk + 1 : nk + nb)].';
      k = k / 2;
      ll = ll + 1;
      if ll == 20
        st = 1;
      end
      if ll == 10
        gndir = Vd / norm(R) * length(R);
        k = 1;
      end
    end
    
    t = t1;
    Vcap = V1;
  end
end
