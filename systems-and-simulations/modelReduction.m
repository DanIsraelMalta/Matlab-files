%{
% modelReduction performs model reduction by removes fast dynamics while maintaining static gain.
% modelReduction is different then matlab sminreal, since it does not look
% for non-minimal states / zero enteries in A,B & C matrix or zero-pole
% cancellation.
%
% sys - state space system
% tol - dynamics faster than this are removed
% sysr - model reduced system
%
% Dan I. Malta 2009
%}
function sysr = modelReduction(sys, tol)
    % extract system characteristics
    [a, b, c, d] = ssdata(sys);
    [U, D] = eig(a);
    
    % housekeeping
    tol = abs(tol);
    n = size(a, 1);
    V = [];
    at = zeros(n);
    
    % run through system matrix
    ii=0;
    while ii < n
        ii = ii + 1;
        
        % real
        if imag(D(ii, ii)) ~= 0
            V = [V real(U(:, ii)) imag(U(:, ii))];
            at(ii, ii) = real(D(ii, ii));
            at(ii + 1, ii) = -imag(D(ii, ii));
            at(ii, ii + 1) = -at(ii + 1, ii);
            at(ii + 1, ii + 1) = at(ii, ii);
             ii = ii + 1;
        else
            V = [V U(:, ii)];
            at(ii, ii) = D(ii, ii);
        end
    end
    
    % rest of matrix
    bt = V \ b;
    ct = c * V;
    dt = d;
    ar = at;
    br = bt;
    cr = ct;
    
    % tolerance check
    jj=0;
    for ii = 1 : n
        jj = jj + 1;
        if abs(at(ii, ii)) > tol
            ar = rowcomp(colcomp(ar, jj), jj);
            br = rowcomp(br, jj);
            cr = colcomp(cr, jj);
            jj = jj - 1;
        end
    end
    
    % D matrix
    G0 = d - c / a * b;
    dr = G0 + cr / ar * br;
    
    % output
    sysr = ss(ar, br, cr, dr);
end

% suppress matrix rows, given their index
function mo = rowcomp(mi, i)
    mo = mi(setdiff(1 : size(mi, 1), i), :);
end

% suppress matrix columns, given their index
function mo = colcomp(mi, i)
    mo = mi(:, setdiff(1 : size(mi, 2), i));
end
