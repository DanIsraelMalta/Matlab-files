%{
% danBuffer partitions a vector into nonoverlapping data segments (frames) of a given length.
% if needed, zero padding is employed.
%
% x - vector
%  n - frame length
% y - reshaped matrix
%
% example:

y = danBuffer(1:15, 6)

% Dan I. Malta 2011
%}
function y = danBuffer(x, n)
    % housekeeping
    len = length(x);
    colmn = ceil(len / n);
    y = zeros(n, colmn);
    
    for i = 1 : colmn
        sInd = max((i - 1) * n + 1, 1);
        eInd = min(sInd + n - 1, len);
        frame = x(sInd : eInd);
        if (1 + eInd - sInd < n)
            frame(end + 1 : end + n - (eInd - sInd) - 1) = 0.0;
        end
        y( :, i) = frame;
    end
end
