%{ 
% envelope retrieve the envelope of a given 2D data set
%
%  x, y - data set
% z - envelope of {x, y} holded in z.top and z.bottom
%
% example:

t = 0 : 0.1 : 4*pi;
y = sin(t) + 0.2 * randn(size(t));
y(end) = 0;
z = envelope(t, y);

figure;
plot(t, y, 'r', t, z.top, 'b', t, z.bottom, 'b');
legend('signal', 'envelope');
grid on;

%
% Dan I. Malta 2008
%}
function z = envelope(x, y)
    % extreme points value & index
    extrMaxValue = y(find(diff(sign(diff(y))) == -2) + 1);
    extrMaxIndex = find(diff(sign(diff(y))) == -2) + 1;
    extrMinValue = y(find(diff(sign(diff(y))) == 2) + 1);
    extrMinIndex = find(diff(sign(diff(y))) == 2) + 1;
    
    % interpolated evnelope
    z.top = pchip(x(extrMaxIndex), extrMaxValue, x);
    z.bottom = pchip(x(extrMinIndex), extrMinValue, x);
end
