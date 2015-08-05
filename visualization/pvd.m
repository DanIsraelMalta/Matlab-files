%{
% pvd performs a progressive vector plot of a given series {x, y}
%
% x, y - data series
% xo, yo - plot initial location
%
% example:

t = 0 : 0.2 : 2 * pi;
x = sin(t);
y = cos(t);
pvd(0, 0, x, y);
hold on;
pvd(6, 4, x, y);
pvd(8, 0, x, y);
pvd(14, 4, x, y);
pvd(16, 0, x, y);
grid on;

%
% Dan I. Malta 2008
%}
function pvd(xo, yo, x, y)
    % housekeeping
    n = length(x);
    innX = isnan(x);
    innY = isnan(y);
    x(innX) = 0;
    y(innY) = 0;
    
    % aggregated  position
    posX = cumsum([xo ; x(:)]);
    posY = cumsum([yo ; y(:)]);
    posX([isnan(1) ; innX(:)]) = NaN;
    posY([isnan(1) ; innY(:)]) = NaN;
    
    % progressive plot
    quiver(posX(1 : n), posY(1 : n), x(:), y(:), 0);
end
