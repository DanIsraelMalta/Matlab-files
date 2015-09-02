%{
% kmeans perform the k-means clustering algorithm using euclidean distance metric
%
% x - a matrix whose rows correspond to points and columns to variables
% k - number of clusters
% tol - tolerance stopage criteria - difference betwen two consecutive runs 'maximal distance' value
%             if not supplied - taken as 1e-6
% maxIter - iteration stopage criteria - maximal number of algorithm iterations
%             if not supplied - taken as 100
% centroids - matrix holding cluster centers
% distLabel - matrix holding the distance from variable point to every
%                            centroid, centroid index it belongs to and
%                            minimum distance among all distance's.
% cause - termination caues (0 - tolerance, 1 - iterations)
%
% example:

% data
l = 2000;
x = [randn(l, 2) + ones(l, 2); randn(l, 2) - ones(l, 2); 0.75 * randn(l, 2) + 3.5 * ones(l, 2)];

% clustering
tic;
k = 3;
[cent, distl, cause] = kmeans(x, k); % 1e-6 tolerance, maximum of 100 iterations
toc;
if cause == 0
    disp('stopage criteria = tolerance.');
else
    disp('stopage criteria = iterations.');
end

% visualization
figure;
cmap = {'b', 'r', 'g', 'm'};
hold on;
for i = 1 : k
    % dots
    plot(x(distl(:, k+1) == i, 1), x(distl(:, k+1) == i, 2), '.', 'color', cmap{i}, 'MarkerSize', 12);

    % closing area
    xx = x(distl(:, k+1) == i, 1);
    yy = x(distl(:, k+1) == i, 2);
    kk = convhull(xx, yy);
    plot(xx(kk), yy(kk), '--', 'color', cmap{i});
end
% centroid's
plot(cent(:, 1), cent(:, 2), 'k*', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'k');
plot(cent(:, 1), cent(:, 2), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
grid on;
title('true centers are @ [-1, -1], [1, 1], [3.5, 3.5]');

%
% Dan I. Malta 2015
%}
function [centroids, distLabel, cause] = kmeans(x, k, tol, maxIter)
    % input argument
    if nargin == 2
        tol = 1e-6;
        maxIter = 100;
    elseif nargin == 3
        maxIter = 100;
    end
        
    % housekeeping
    cause = 1;
    len = size(x, 1);
    distLabel = zeros(size(x, 1), k + 2);
    maxDistancePrev = 0;

    % randomaly (uniformaly) place centroids along the data
    centroids = x(ceil(rand(k, 1) * len), :);
    
    % iterate away
    for iter = 1 : maxIter
        % distance'ing & labeling
        for i = 1 : len
            % euclidean distance calculation
            for j = 1 : k
                distLabel(i, j) = norm(x(i, :) - centroids(j, :));
            end

            % minimal distance and cluster index
            [distLabel(i, k + 2), distLabel(i, k + 1)] = min(distLabel(i, 1 : k));
        end
        
         % tolerance stopage criteria
         [maxDistance, dummyCenter] = max(distLabel(i, 1 : k));
         if abs(maxDistance - maxDistancePrev) <= tol
             cause = 0;
             return;
         end
         maxDistancePrev = maxDistance;
        
        % clustring & cenroid'ing
        for i = 1 : k
            % k clusters
            A = (distLabel(:, k + 1) == i);
            
            % new cluster centers
            centroids(i, :) = mean(x(A, :));
        end
    end
end
