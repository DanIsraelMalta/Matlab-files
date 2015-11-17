%{
% logSpan logarithmically grid's a given region with required number of points
%
% vmin - gridded region minimal value [vmin, ..]
% vmax - gridded region maximal value [..., vmax]
% n - grid size
% span - [vmin, vmax] logarithmically gridded with n points
%
% Dan I. Malta 20115
%}
function span = logSpan (vmin, vmax, n)
    fac = exp((log(vmax) - log(vmin)) / (n - 1));
    w = ones(1, n) * fac;
    span = cumprod(w) * vmin / fac;
end
