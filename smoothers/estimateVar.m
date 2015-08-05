%{
% estimateVar estimates the variance of an additive noise upon a given signel.
% this is done by "fitting" a smooth "thin-plate" spline thru the data and
% then using a maximum apriori variance estimator.
%
% unlike the sampled variance calculation, which maximizes the maxmimum
% likelihood, this algorithm  maximizes the aposteriory probability,
% thus it is locally and asymptotically unbiased! 
%
% y - observation series
% noisevar - estimated noise variance
%
% example:

Ts=0.05;
time = linspace(0,8*pi,Ts*2500);
y = cos(time/pi)+(time.^0.57)/3 - 1.0;
var0 = 0.03; % noise variance
yn = y + sqrt(var0)*randn(size(y));
v=estimateVar(yn); % noise estimated variance

% filter the noise accroding to the variance (first order FIR, noise variance is filter time constant)
yf = filter([Ts Ts], [Ts+2*v Ts-2*v], y);

% plot
figure('color',[1 1 1]);
plot(time,yn,'g',time,y,'b',time,yf,'r'); grid on;
xlabel('Time [sec]');ylabel('Amplitude');
legend('noisy signal [cos(t/\pi)+0.3*t^{0.57}-1+N(0, 0.03)]', 'original signal [cos(t/\pi)+0.3*t^{0.57}-1]', 'filtered noisy signal');
title(['Noise variance is estimated as ', num2str(v),...
                   ' (true noise variance is ', num2str(var0), '); RMS difference betwen filtered signal and original signal is ', num2str(norm(y-yf)/sqrt(1/Ts))]);

%
% Dan I. Malta 2010
%}
function noisevar = estimateVar(y)
    % housekeeping
    d = length(size(y));
    siz = size(y);
    S = zeros(siz);
    
    % "thin plate" spline
    for i = 1:d
        siz0 = ones(1,d);
        siz0(i) = siz(i);
        S = bsxfun(@plus,S,cos(pi*(reshape(1:siz(i),siz0)-1)/siz(i)));
    end
    S = 2*(d-S(:));

    % N-D Discrete Cosine Transform of Y
    y = dctn(y);
    y = y(:);

    % square calculations (better performance)
    S = S.^2; y = y.^2;

    % Upper and lower bounds for the smoothness parameter
    N = sum(siz~=1); % tensor rank of the y-array
    hMin = 1e-6; hMax = 0.99;
    sMinBnd = (((1+sqrt(1+8*hMax.^(2/N)))/4./hMax.^(2/N)).^2-1)/16;
    sMaxBnd = (((1+sqrt(1+8*hMin.^(2/N)))/4./hMin.^(2/N)).^2-1)/16;

    % Minimization of the GCV score
    fminbnd(@func,log10(sMinBnd),log10(sMaxBnd),optimset('TolX',.1));

    % Generalized cross validation score
    function score = func(L)
        M = 1-1./(1+10^L*S);
        noisevar = mean(y.*M.^2);
        score = noisevar/mean(M)^2;
    end

    % discrete cosine transform
    function y = dctn(y)
        y = double(y);
        sizy = size(y);
        y = squeeze(y);

        dimy = ndims(y);
        if isvector(y)
            dimy = 1;
            y = y(:).';
        end

        if ~isreal(y)
            y = complex(dctn(real(y)),dctn(imag(y)));
        else
            for dim = 1:dimy
                y = shiftdim(y,1);
                siz = size(y);
                n = siz(1);

                y = y([1:2:n 2*floor(n/2):-2:2],:);
                y = reshape(y,n,[]);

                w = exp(1i*(0:n-1)'*pi/2/n);
                y = y*sqrt(2*n);
                y = ifft(y,[],1);
                y = bsxfun(@times,y,w);
                y = real(y);
                y(1,:) = y(1,:)/sqrt(2);

                y = reshape(y,siz);
            end
        end

        y = reshape(y,sizy);
    end
end
