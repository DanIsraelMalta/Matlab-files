%{
% danDct performs the direct or inverse discrete cosine transform on a given vector.
%
% x - vector to be transformed
% type - a string defining the type of conversion. options are: 'direct' / 'inverse'
% y - cosine transform output
%
% Dan I. Malta 2013
%}
function y = danDct(x, type)
    % housekeeping
    len = numel(x);
    x = x(:);
    
    % conversion
    switch type
        case 'inverse'
            % weights
            w = len * exp(i * (0 : len - 1) * pi / (2 * len)).';
            
            % cosine transformation
            ydct = real(ifft(w .* x));
            
            % output
            y = zeros(size(x));
            y(1 : 2 : end) = ydct(1 : len / 2);
            y(2 : 2 : len) = ydct(len : -1 : len / 2 + 1);
        case 'direct'
              % weights
              w = [1; 2 * (exp(-i * (1 : len - 1) * pi / (2 * len))).'];
              
              % input re-ordering
              x = [x(1 : 2 : end, :); x(end : -2 : 2, :)];
              
              % output
              y = real(w .* fft(x));
        otherwise
            y = 0;
    end
end
