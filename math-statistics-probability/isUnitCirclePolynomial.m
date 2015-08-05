%{
% isUnitCirclePolynomial test if a given polynomial roots lie inside the unit circle.
% 
% This solution is better then matlab's "all(abs(roots(C))<1)" command,
% because it is much-much faster (4-60 times faster, since it doesn't
% compute the companion matrix eigenvalues).
%
% c - coefficients of polynomial c(1)*x^n + ... + c(n)*x + c(n+1)
% b - '1' if polynom roots are within the unit circle, '0' - otherwise
 %
% example:
%

numOfRoots = 500;
pnotstable = poly(1.1 * rand(1, numOfRoots));  % unit circle not stable polynomial
pstable = poly(rand(1, numOfRoots) / numOfRoots);  % unit circle stable polynomial

clc;
disp(['evaluation of two polynoms with ', num2str(numOfRoots), ' roots (first one is not stable, second one is stable).']);
disp('---');

disp('isUnitCirclePolynomial evaluation:');

tic;
b = isUnitCirclePolynomial(pnotstable);
t = toc;
disp(['polynom is ', num2str(b), ' (1 = stable, 0 = not stable), evalutaion time: ', num2str(t)]);

tic;
b = isUnitCirclePolynomial(pstable);
t = toc;
disp(['polynom is ', num2str(b), ' (1 = stable, 0 = not stable), evalutaion time: ', num2str(t)]);

disp('---');
disp('matlab evaluation:');
tic;
b = all(abs(roots(pnotstable))<1);
t = toc;
disp(['polynom is ', num2str(b), ' (1 = stable, 0 = not stable), evalutaion time: ', num2str(t)]);

tic;
b = all(abs(roots(pstable))<1);
t = toc;
disp(['polynom is ', num2str(b), ' (1 = stable, 0 = not stable), evalutaion time: ', num2str(t)]);

%
% Dan I. Malta 2012
%}
function b = isUnitCirclePolynomial(c)
    % housekeeping
    [lr,lc] = size(c);
    
    % JURY-Scheme
    b = ones(lr, 1);
    lambda = zeros(lr, 1);
    while lc > 1
        lambda = c(:, lc) ./ c(:,1);
        b = b & (abs(lambda) < 1);
        
        % reduced polynomial must be stable as well
        c(:, 1 : lc - 1) = c(:, 1 : lc - 1) - lambda(:, ones(1, lc - 1)) .* c(:, lc : -1 : 2);
        lc = lc - 1;
    end
end

