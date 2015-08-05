%{
% kolmogorovSmirnov performs the two-sample kolmogorov-smirnov goodnes
% of fit hypothesis test - determine in two independent random samples are
% drawn from the same underlying population.
%
% For years I mocked the need of this test, until I read Yakov Bar Shalom's
% "Estimation with applications to tracking and navigation".
% I built this function to have the same API as in kstest2 from the
% statistical toolbox, which off course - I do not have.
%
% x1, x2 - samples
% alpha - desired significance level (choose 0.05 for 95%)
% tail - test type (default should be 'unequal')
%        'unequal' - two sided test
%        'larger'  - one sided test (cdf(x1) > cdf(x2))
%        'smaller' - one sided test (cdf(x1) < cdf(x2))
% H - hypothesis test
%     H = 0 -> do not reject thr null hypothesis at significande level alpha
%     H = 1 -> reject the null hypothesis at significande level alpha
% p - asymptotic p value (alpha > p -> reject hypothesis)
% ks - max|empirical distribution function of x1 - empirical distribution function of x2|
%
% Dan I. Malta 2010
%}
function [H, pValue, KSstatistic] = kolmogorovSmirnov(x1, x2, alpha, tail)
  % columnize and NaN'ize
  x1 = x1(:);
  x2 = x2(:);
  x1 = x1(~isnan(x1));
  x2 = x2(~isnan(x2));

  % estimate empirical distribution function
  binE = [-inf ; sort([x1; x2]); inf];
  binC1 = histc(x1, binE, 1);
  binC2 = histc(x2, binE, 1);
  sumC1 = cumsum(binC1) ./ sum(binC1);
  sumC2 = cumsum(binC2) ./ sum(binC2);
  sCDF1 = sumC1(1 : end - 1);
  sCDF2 = sumC2(1 : end - 1);
  
  % kolmogorov-smirnov statistics
  switch tail
    case 'unequal'
      dCDF = abs(sCDF1 - sCDF2);
    case 'larger'
      dCDF = sCDF2 - sCDF1;
    case 'smaller'
      dCDF  =  sCDF1 - sCDF2;
  end
  KSstatistic = max(dCDF);
  
  % asymptotic P-value
  n1 = length(x1);
  n2 = length(x2);
  n = n1 * n2 / (n1 + n2);
  lambda = max((sqrt(n) + 0.12 + 0.11 / sqrt(n)) * KSstatistic , 0);
  if ~strcmp(tail, 'unequal')
    pValue = exp(-2 * lambda * lambda);
  else
    j =  (1 : 101)';
    pValue = 2 * sum((-1) .^ (j - 1) .* exp(-2 * lambda * lambda * j.^2));
    pValue = min(max(pValue, 0), 1);
  end
  
  % hypothesis
  H = alpha >= pValue;
end
