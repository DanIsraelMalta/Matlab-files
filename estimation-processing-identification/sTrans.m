%{
% sTrans implements the Stockwell transform (with out the use of loops)
%
% x - signal
% st - stockwell transform of x
%
% Dan I. Malta 2011
%}
function st = sTrans(x)
    % housekeeping
    [m, N] = size(x); %#ok<ASGLU>
    nxaf = fix(N / 2);
    odvn=1;
    
    % half the length?
    if 2 * nxaf == N
        odvn=0;
    end
    
    % fft
    f = [0 : nxaf -nxaf + 1 - odvn : -1 ]/ N;
    xft=fft(x);

    % calculate the gaussians as one matrix
    invfk = 1 ./ f(2 : nxaf + 1);
    invfk = invfk';
    W = 2 * pi * repmat(f, nxaf, 1) .* repmat(invfk, 1, N);
    G=exp((-W.^2) / 2);

    % calculate the toeplitz matrix
    xW = toeplitz(xft(1 : nxaf + 1)', xft);
    
    % remove zero frequency (first row)
    xW = xW(2 : nxaf + 1, :);

    
    % stockwell Transform
    st = ifft(xW .* G, [], 2);

    %Add zero frequency
    st0 = mean(x) * ones(1,N);
    st=[st0;st];
end
