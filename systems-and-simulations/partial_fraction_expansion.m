%{
% partial_fraction_expansion finds the residues, poles and direct term of a
% partial fraction expansion of a given set of zeros, poles and a gain.
%
%       Z(s)      R(1)       R(2)            R(n)
% Kin * ---- => -------- + -------- + ... + -------- + k(s)
%       P(s)    s - p(1)   s - p(2)         s - p(n) 
%
% input: Kin, Z, P
% output: R, p, k
%
% Remarks:
% > k is empty if length(Z) < length(P), otherwise, for a proper system: length(k) = 1
% > if P(j) = ... = P(j + m - 1) is a pole of multiplicity m, then the
%    expansion includes terms of the form:
%                R(j)        R(j+1)                R(j+m-1)
%              -------- + ------------   + ... + ------------
%              s - p(j)   (s - p(j))^2           (s - p(j))^m
% > the function does not work with repeated complex poles!
%
% Dan I. Malta (2016)
%}
function [xo_r, xo_p, xo_k] = partial_fraction_expansion(xi_z, xi_p, xi_kin1)
    % housekeeping
    xo_r = zeros(length(xi_p), 1);
    xo_k = [];
    if length(xi_z) == length(xi_p)
        xo_k = 1;
    end
    
    % Q matrix (poles and their conjugates)
    [~, i1] = sort(-abs(xi_p));
    xi_p    = xi_p(i1);
    pLength = length(xi_p);
    Q       = zeros(1,pLength);
    
    ct = 1;
    while ct <= pLength
        rt       = sort(find(xi_p(ct) == xi_p));
        rtLength = length(rt);
        Q(rt)    = (rtLength-1) : -1 : 0;
        
        if all(diff(rt) == 0)
            ct = ct + rtLength;
        else
            ct = ct + 1;
        end
    end
    Q = [conj(xi_p(:)'); Q];
 
    % residue
    for j1 = 1 : length(xi_p)
        p1 = Q(1,j1);
        n  = Q(2,j1);
        temp = xi_p(xi_p ~= p1);
        
        if ~isempty(xo_k)
            tmp = xo_k * prod(-xi_p + p1);
        else
            tmp = 0;
        end
        
        if n == 0
            xo_r(j1) = (prod(-xi_z + p1) - tmp) / prod(-temp + p1);
        else
            v = poly(temp);
            u = poly(xi_z) - tmp;
            
            for j = 1 : n
                [u,v] = polyder(u, v);
            end
            
            c = 1;
            if xo_k < n
                c = prod(1 : n - xo_k);
            end
            
            xo_r(j1) = polyval(u, p1) / prod(-temp + p1)^(2 ^ n) / c;
        end
    end
 
    % poles
    xo_p = xi_p(:);
    
    % scale by input gain
    xo_k = xi_kin1 * xo_k;
    xo_r = xi_kin1 * xo_r;
    
    % are complex pairs also conjugate?
    indi1 = find(imag(xi_p) > 0);
    indi2 = find(imag(xi_p) < 0);
    if ~isempty(indi1)
        for j = 1 : length(indi1),
            xo_r(indi2(j)) = xo_r(indi1(j))';
        end
    end
    
    % make sure residues of real poles are zero
    ind = find(imag(xi_p) == 0);
    if ~isempty(ind)
        xo_r(ind) = real(xo_r(ind));
    end
end
