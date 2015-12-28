%{
% space(start, end, len, type) generates a vector with specially spaced 'len'
% points in the range [start, end], according to the 'type' parameter (a string):
%  type                          spaced vector
% ------                       -----------------
% "EDGE"                      start    || |  |   |    |     |     |    |   |  | ||  end
% "CENTER"                    start |     |    |   |  | | | | | |  |   |    |     | end
% "END"                       start |      |      |      |     |    |   |  | ||     end
% "START"                     start || |  |   |    |     |      |      |      |     end
%
%
% default are 100 points & 'EDGE' type spacing.
%
% example:

len = 60;
yEdge   = space(0, 100, len, 'EDGE');
yCenter = space(0, 100, len, 'CENTER');
yEnd    = space(0, 100, len, 'END');
yStart  = space(0, 100, len, 'START');

figure;
subplot(2,2,1);
hist(yEdge, len);
title('edge');
subplot(2,2,2);
hist(yCenter, len);
title('center');
subplot(2,2,3);
hist(yEnd, len);
title('end');
subplot(2,2,4);
hist(yStart, len);
title('start');

% Dan I. Malta 2015
%}
function xo_vec = space(xi_start, xi_end, xi_length, xi_type)
    % housekeeping
    if nargin == 2
        xi_length = 100;
        xi_type = 'EDGE';
    elseif nargin == 3
        xi_type = 'EDGE';
    end
    xo_vec = zeros(xi_length, 1);
    
    % create spaced vector
    switch xi_type
        case 'EDGE'
            xi_length = floor(double(xi_length));
            xi_length = xi_length - 1;
            xo_vec = xi_start + (1 - cos(pi / xi_length * (0 : xi_length))) * (xi_end - xi_start) / 2;
        case 'CENTER'
            xi_length = floor(double(-xi_length));
            xi_length = xi_length + 1;
            xo_vec = xi_start + acos(1 - 2 * (0 : -1 : xi_length) / xi_length) * (xi_end - xi_start) / pi;
        case 'END'
            xi_length = floor(double(xi_length));
            xi_length = xi_length - 1;
            xo_vec = xi_start + sin(pi / (2 * xi_length) * (0 : xi_length)) * (xi_end - xi_start);
        case 'START'
            xi_length = floor(double(-xi_length));
            xi_length = xi_length + 1;
            xo_vec = xi_end + sin(pi / 2 * (1 - (0 : -1 : xi_length) / xi_length)) * (xi_start - xi_end);
        otherwise
            error('space.m: unknown spacing type.'); %#ok<ERTAG>
    end
end
