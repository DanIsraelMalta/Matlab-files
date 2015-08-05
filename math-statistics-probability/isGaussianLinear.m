%{
% isGaussianLinear tests a time series for gaussianity (no skewness) and linearity.
%
% x - time series
% c - resolution (smoothing) parameter (0.5, 1.0);
%       increasing c reduces the variance but increases the bias (window bandwidth)
% nfft - fft length
% sg - statistics of gaussianity test (1x3 vector(
%         [observed chi square value,
%         degrees of freedom (of chi square value,
%         probability that test is wrong]
% sl - statistics of linearity test (1x3 vector(
%         [estimated inter-quantile range,
%         scaled version of mean skewness,
%         theoretical inter-quantile range]
%
% Dan I. Malta 2015
%}
function [sg, sl] = isGaussianLinear(x, c, nfft)
    % housekeeping
    nsamp = length(x);
    if nfft  > nsamp
        error('fft length is longer then sample length.') %#ok<ERTAG>
    end
    nrecs = floor(nsamp / nfft);
    nrecs = max(nrecs, 1);
    ksamp = min(nfft, nsamp);
    c = min(1, max(0.5, c));
    M = fix(nfft ^ c);
    M = M + 1 - rem(M, 2);
    Msmuth = M;
    
    % raw bi-spectrum (F) and power spectrum (S) estimation
    mrow = fix(nfft / 3) + 1;
    ncol = fix(nfft / 2) + 1;
    F = zeros(mrow, ncol);
    S = zeros(nfft, 1);
    mask  = hankel([1 : mrow], [mrow : mrow + ncol - 1]); %#ok<NBRAK>
    ind = 1 : ksamp;
    for k = 1 : nrecs
        y = x(ind);
        xf = fft(y - mean(y), nfft);
        xc = conj(xf);
        
        % power spectrum
        S = S + xf .* xc;
        F = F + xf(1 : mrow) * xf(1:ncol).' .* reshape(xc(mask), mrow, ncol) ;
        ind = ind + ksamp;
    end
    F = F / (nfft * nrecs);
    S = S / (nfft * nrecs);
    
    % nullify area outsie principle domain
    ind = (0 : (mrow - 1))';
    F(1 : mrow, 1 : mrow) = triu(F(1 : mrow, 1 : mrow));
    Q = ones(mrow, ncol);
    Q(1 : mrow, 1 : mrow) = triu(Q(1 : mrow, 1 : mrow)) + eye(mrow);
    r = (rem(nfft, 3) == 2);
    for k = mrow + 1 : ncol - r
        j = k - mrow;
        Q(mrow - 2 * j + 2 : mrow, k) = zeros(2 * j - 1, 1);
        F(:, k) = F(:, k) .* Q(:, k);
        Q(mrow - 2 * j + 1, k) = 2;
    end
    
    % remove DC term
    F = F(2 : mrow, 2 : ncol);
    Q = Q(2 : mrow, 2 : ncol);
    mrow = mrow - 1;
    ncol = ncol - 1;
    
    % smooth estimated power spectrum / bispectrum
    M = Msmuth;
    mby2 = (M + 1) / 2;
    m1 = rem(mrow, M);
    m2 = rem(ncol, M);
    m1 = - m1 + M * (m1 ~= 0);
    m2 = - m2 + M * (m2 ~= 0);
    F = [F, zeros(mrow, m2); zeros(m1, ncol + m2)];
    Q = [Q, zeros(mrow, m2); zeros(m1, ncol + m2)];
    k  = size(F) / Msmuth;
    k1 = k(1);
    k2 = k(2);
    
    % kronecker tensor
    ind = 1 : Msmuth;
    B = zeros(k1, k2);
    Q1 = B;
    for i = 1 : k1
        for j = 1 : k2
            t = F( (i - 1) * Msmuth + ind, (j - 1) * Msmuth + ind);
            B(i, j) = mean(t(:));
            t = Q( (i - 1) * Msmuth + ind, (j - 1) * Msmuth + ind);
            Q1(i, j) = sum(t(:));
        end
    end
    Q = Q1;
    
    % box-car smoother
    M = Msmuth;
    S = conv(S, ones(M, 1)) / M;
    S = S(M + 1 : M : M + nfft - 1);
    S1 = B*0;
    S1 = S(1 : k1) * S(1 : k2)' .* hankel(S(2 : k1 + 1), S(k1 + 1 : k1 + k2));
    S1 = ones(k1, k2) ./ S1;
    ind = find(Q > 0);
    Q   = Q(ind);
    Bic = S1(ind) .* abs(B(ind)).^2;
    
    % gaussianity (non-skewness) test
    scale = 2 * (Msmuth^4) / nfft;
    Xstat = scale * Bic ./ Q;
    df = length(Q) * 2;                          % degrees of freedom
    chi_val = sum(Xstat);                % observed value of chi-square r.v.

    % probability that gaussianity test is wrong
    h = 1/3;
    b = sqrt(2 / df) / 3;
    a = (1 - b^2);
    y = (chi_val / df)^(1/3);
    y = (y - a) / b;
    pfa = 0.5 * erfc(y / sqrt(2));
    pfa = round(pfa * 10000) / 10000;     % five decimal places
    
    sg = [chi_val, df, pfa];
    
    % linearity test
    
    % first and third quantile search
    ind = find(Q == M^2);
    rtest = Xstat(ind); %#ok<FNDSB>
    rtest = sort(rtest);
    l1 = length(rtest);
    if l1 < 4
        error('time series is to long or fft length to short for linearity test'); %#ok<ERTAG>
    else
        lo1 = fix(l1 / 4);
        lo2 = fix(l1 / 4+0.8);
         hi1 = fix(3 * l1 / 4);
         hi2 = fix(3 * l1 / 4+0.8);
         r1  = (rtest(lo1) + rtest(lo2)) / 2;
         r3  = (rtest(hi1) + rtest(hi2)) / 2;
         Rhat  = r3 - r1;
    end
    
    % estimate mean skewness
    lam = mean(Bic) * 2 * Msmuth^2 / nfft;
    
    % theoretical inter-quantile range of chi2(f=2, lam)
    df = 2;
    h  = 1/3 + (2/3) * ( lam / (df + 2 *lam) )^2;
    a = 1 + h * (h - 1) * (df + 2*lam) / ( (df + lam)^2) - h * (h-1) * (2-h) * (1-3*h) * (df + 2*lam)^2 / (2 * (df+lam)^4);
    b = h * sqrt(2 * (df+2*lam) ) / (df+lam)  * (  1 - (1-h) * (1-3*h) * (df+2*lam) / ( 2*(df+lam)^2) );
    yn   = sqrt(2) * erfinv(2*[-0.25,.25]);       % quartiles of Phi
    xc   = (a + b * yn) .^ (1/h) * (df + lam);         % quartiles of Chi
    Rth  = xc(2)-xc(1);                          % IQ range
    
    sl = [Rhat, lam, Rth];
end
