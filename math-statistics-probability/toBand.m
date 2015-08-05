%{
% toBand expands a given matrix to a banded matrix of a given size, where
% the input matrix is it's local operator.
%
% D - input matrix
% n - size of banded matrix (nxn)
% s - banded matrix (nxn) where D acts as the local operator
%
% example:

S = toBand(3 * ones(5) - 2 * eye(5), 10)

% Dan I. Matlab 2013
%}
function S = toBand( D, n)
    % housekeeping
    [p, q] = size(D);
    
    % Is S sufficiently large
    if n <= max( [p,q])
        error( 'n > max( [p,q] ).'); %#ok<ERTAG>
    end
    
    % if q is even then the computation is valid at the interstitial points (consequently p should be odd)
    if ~rem(q, 2)
        if p ~= (q - 1)
            error('Given [p,q] = size( D ), if q is even then p = q - 1.'); %#ok<ERTAG>
        end
        k = round((p - 1) / 2);
        m = n - 1;
    else  % if q is odd then the computation is valid at the nodes
        if p ~= q
            error('Given [p,q] = size( D ), if q is odd then p = q.'); %#ok<ERTAG>
        end
        k = round((p - 1) /2);
        m = n;
    end

    % Extract the partitioning of D
    D1 = D(1 : k, :);
    d = D(k + 1, :)';
    D2 = D(k + 2 : end, :);
    rd = 1 : length(d);
    
    % housekeeping (2)
    S = zeros( m, n);
    
    % add the start and end blocks
    rs = 1 : k;
    cs = 1 : q;
    [Cs, Rs] = meshgrid(cs, rs);
    R = Rs(:);
    C = Cs(:);
    D = D1(:);
    rs = (m - k + 1) : m;
    cs = (n - q + 1) : n;
    [Cs, Rs] = meshgrid(cs, rs);
    R = [R; Rs(:)];
    C = [C; Cs(:)];
    D = [D; D2(:)];
    
    % add diagonal vector
    noVs = m - 2 * k;
    for j = 1 : noVs
        rs = j + k;
        cs = j + rd - 1;
        [Rs, Cs] = meshgrid(rs, cs);
        R = [R; Rs(:)];
        C = [C; Cs(:)];
        D = [D; d];
    end
    
    % generate S
    S = sparse(R, C, D, m, n);
    
    % full matrinx conversion
    S = full( S ); %#ok<ACCUM>
end
