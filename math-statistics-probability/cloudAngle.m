%{
% cloudAngle estimate the main axis angle of a point cloud using principle component analysis
%
% x - data vector
% y - data vector
% theta - main axis angle
%
% example:

np = 3000;
angle = 38 * pi / 180;
dcm = [cos(-angle), sin(-angle);
              -sin(angle), cos(-angle)];
x = 4 * randn(1, np);
y = 2*randn(1, np);
rotxy = dcm * [x; y];
alpha = cloudAngle(rotxy(1, :), rotxy(2, :));

figure;
plot(rotxy(1, :), rotxy(2, :), '.r')
axis equal;
hold on;
quiver(0, 0, cos(alpha) * norm(axis), sin(alpha) * norm(axis), 'k', 'linewidth', 2);
title(['cloud angle is ', num2str(angle * 180 / pi), ', estimated angle is ', num2str(180 - alpha * 180 / pi)]);
grid on;

%
% Dan I. Malta 2010
%}
function theta = cloudAngle(x, y)
    % detrend and columnize
    x = x(:) - mean(x(:));
    y = y(:) - mean(y(:));
    
    % principle component analysis
    c = cov(x, y);
    [a, ev] = eig(c);
    [ev, ind] = sort(diag(ev)); %#ok<ASGLU>
    
    % main axis angle
    [xa, ya] = deal(a(1, ind(end)), a(2, ind(end)));
    theta = cart2pol(xa, ya);
end
