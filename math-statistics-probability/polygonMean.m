%{
% polygonMean accepts any number of 2D polygons (given on a cartesian canvas)
% and return their average (mean) polygonial shape, centroid and area.
%
%
% remarks:
% 1) after a lot of hours using this function, I can say that since polygonMean assumes
%    "not to iregular" polygons, which means that polygons with lines radiating outward
%    from ceteroid (or something like that) will not be traeted correctly.
%    otherwise - the function works very well.
% 2) the algorithm I use perform polar-cartesian transformation (forward & backward)
%    and that makes it sensitive to rounding errors.
%
% [x, y, cx, cy, a] = polygonMean(x1, y1, x2, y2, ...)
%
% example:

% angles
th = linspace(0, 2*pi, 100);

figure;

% four covairance ellipse
x1=0; y1=0; tilt= 30 * pi / 180; rx = 10; ry = 20;
plot(x1, y1, '*k');
hold on;
x1 = x1 + rx*cos(th)*cos(tilt) - ry*sin(th)*sin(tilt);
y1 = y1 + rx*cos(th)*sin(tilt) + ry*sin(th)*cos(tilt);

x2=0.5; y2=0.5; tilt= 35 * pi / 180; rx = 12; ry = 18;
plot(x2, y2, '*k');
x2 = x2 + rx*cos(th)*cos(tilt) - ry*sin(th)*sin(tilt);
y2 = y2 + rx*cos(th)*sin(tilt) + ry*sin(th)*cos(tilt);

x3=2; y3=2; tilt= 40 * pi / 180; rx = 8; ry = 20;
plot(x3, y3, '*k');
x3 = x3 + rx*cos(th)*cos(tilt) - ry*sin(th)*sin(tilt);
y3 = y3 + rx*cos(th)*sin(tilt) + ry*sin(th)*cos(tilt);

x4=2; y4=0.5; tilt= 25 * pi / 180; rx = 10; ry = 18;
plot(x4, y4, '*k');
x4 = x4 + rx*cos(th)*cos(tilt) - ry*sin(th)*sin(tilt);
y4 = y4 + rx*cos(th)*sin(tilt) + ry*sin(th)*cos(tilt);

% average covariance ellipse
[xm, ym, cx, cy, a]=polygonMean(x1', y1', x2', y2', x3', y3');

% mean covariance ellipse angle
alpha = cloudAngle(xm, ym);

plot(x1, y1, 'k', x2, y2, 'k', x3, y3, 'k', cx, cy, '*b');
hold on;
plot(xm, ym, 'b', 'LineWidth', 2);
quiver(cx, cy, cos(alpha) * norm(axis), sin(alpha) * norm(axis), 'b', 'linewidth', 2);
title(['mean polygon area is ', num2str(a), ', its center is at (', num2str(cx), ', ', num2str(cy), '), and its main direction is ', num2str(360 - alpha * 180 / pi), '\circ']);
grid on;

%
% Dan I. Malta 2013
%}
function [x, y, cx, cy, a]=polygonMean(varargin)
  % housekeeping
  centroids = [];
  areas = [];
  radii = zeros(360, 1);
  quantity = length(varargin) ./ 2; % number of polygons (seperated to x, y vectors - must be an even number)

  % iterate over all polygons
  for n = 1 : quantity
    % extract polygon series
    x = varargin{2 * n - 1};
    y = varargin{2 * n};
    
    %calculate centroid
    a = 0;
    cx = 0;
    cy = 0;
    len = length(x);
    for i = 1 : len-1
        common = x(i) * y(i + 1) - x(i + 1) * y(i);
        a = a + common;
        cx = cx + (x(i) + x(i + 1)) * common;
        cy = cy + (y(i) + y(i + 1)) * common;
    end
    a = 0.5 * a;
    cx = cx / (6 * a);
    cy = cy / (6 * a);
    centroids = [centroids; cx cy];
    areas = [areas, a];

    % "detrend" centroid to its origin
    x = x - cx;
    y = y - cy;
    
    % convert to polar coordinates
    [theta, r] = cart2pol(x, y);
    
    % wrap to [0 360]
    theta = round(theta * (180 / pi));
    theta = theta + 180;
    
    % remove zero angles
    a = find(theta == 0);
    theta(a) = [];
    r(a) = [];
    
    % sort angles in ascending order
    [theta, i] = sort(theta);
    r=r(i);
    
    % remove consecutive multiple instances of the same angle
    d = diff(theta);
    a = find(d == 0);
    theta(a) = [];
    r(a) = [];

    % extrapolate (cubic wise) so angle range would apan 1 to 360
    if range(theta) < 359
      while theta(1) ~= 1
        yi = interp1(theta, r, (theta(1) - 1), 'cubic', 'extrap'); 
        theta = [theta(1) - 1 ; theta];
        r = [yi; r];
      end
      
      while theta(end) ~= 360
        yi = interp1(theta, r, (theta(end) + 1), 'cubic', 'extrap'); 
        theta = [theta; theta(end) + 1];
        r=[r; yi];
      end
    end

    % interpolate (linearly) missing angles
    n = 1;
    while length(theta) <= 359
      if (theta(n + 1) - theta(n)) > 1
        yi = interp1(theta, r, n+1);
        theta = [theta(1 : n); n + 1; theta(n + 1 : end)];
        r = [r(1 : n); yi; r(n + 1 : end)];
      end
      n = n + 1;
    end
    radii = [radii r];
  end

  % remove unnecessary zeros
  radii(:, 1) = [];
  
  % mean radii
  radii = mean(radii'); %#ok<UDIM>
  angles = ((1 : 360) - 180) / (180 / pi);
  
  % convert back to cartesian coordinates.
  [x, y] = pol2cart(angles, radii);
  
  % mean area centroid
  a = mean(areas);
  centroids = mean(centroids);
  cx = centroids(1);
  cy = centroids(2);
  
  % trend mean polygon
  x = x + cx;
  y = y + cy;
end
