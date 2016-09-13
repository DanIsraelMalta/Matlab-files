%{
% stepFinder extracts a piecewise constant signal from a noisy observation
% using a maximum likelihood based penalization schema.
%
% notice - this algorithm executes an ietartive solution, therefor it is scaled according to input length.
% so, for high volume data, I strongly suggest to start with a small number of iteration and
% increase that number according to the result.
%
% xi_y        - observation vector from which the piecewise constant signal is extracted
% xi_gamma    - regularization parameter (positive)
% xi_max_iter - maximum number of iterations (defaults is 50)
% xo_y        - piecewise constant signal
%
% example:

signal = [5*ones(1, 50), 7.5 * ones(1, 100), 10 * ones(1, 50), 2.5 * ones(1,150), zeros(1, 100), -3 * ones(1, 100)];
observation = signal + randn(size(signal));
tic; yf = stepFinder(observation, 1.0, 10); time = toc;

figure;
plot(observation, '-.k');
hold all;
plot(signal, 'r');
plot(yf, 'b');
grid on;
legend('observaton', 'signal', 'signal estimate');
title(['underlying signal estimation took ', num2str(time), ' [sec]']);

%
% Dan I. Malta (2016)
%}
function xo_y = stepFinder(xi_y, xi_gamma, xi_max_iter)
    % input arguments
    if (nargin < 3)
        xi_max_iter = 50;
    end
    if (xi_gamma < 0)
        error('stepFinder: regularization parameter must be positive'); %#ok<ERTAG>
    end
    
    % housekeeping
    xi_y = xi_y(:);
    len = length(xi_y);
    Knot = zeros(0,1);      % piecewise knot location
    iter = 1;
    
    % iteratve solution
    while (iter < xi_max_iter)
        % locate new knots
        newKnots = zeros(len, 1);
        for i = 1 : len
            % new knot
            knotI = Knot;
            knotI(end + 1) = i;
            knotI = sort(knotI);

            % optimum knot location
            xtest = zeros(len, 1);
            l = 1 : (knotI(1) - 1);
            xtest(l) = median(xi_y(l));
            for j = 2:iter
                l = knotI(j - 1) : (knotI(j) - 1);
                xtest(l) = median(xi_y(l));
            end
            l = knotI(iter) : len;
            xtest(l) = median(xi_y(l));

            % likelihood of knot location
            newKnots(i) = sum(abs(xtest - xi_y));
        end
        
        % Choose new location that minimizes the likelihood term
        [v, i] = min(newKnots); %#ok<ASGLU>

        % Add to knot locations
        Knot(end + 1) = i;
        Knot = sort(Knot);

        % recalcualte last iteration based upon new data
        i = 1 : (Knot(1) - 1);
        xo_y = zeros(len, 1);
        xo_y(i) = median(xi_y(i));
        for j = 2 : iter
            i = Knot(j - 1) : (Knot(j) - 1);
            xo_y(i) = median(xi_y(i));
        end
        i = Knot(iter) : len;
        xo_y(i) = median(xi_y(i));

        % iteration update
        iter = iter + 1;
    end
end
