%{
% matSymmFac accepts a square matrix and performs the following symmetric factorization:
% P = PP * PP' - PN * PN'
%
% Dan I. Malta (2017)
%}
function [xo_pp, xo_pn] = matSymmFac(xi_p)
    % housekeeping
    [nr, nc]=size(xi_p);
    if nr~=nc
        error('Input matrix not square.');
    else
        n = nr;
    end
    
    % singular value decompsition (P = U*S*V' = u1*s1*v1' + u2*s2*v2' + ...)
    [U, S, V]=svd(xi_p);
    
    % zeros on diagonal
    s=diag(S);
    if0 = find(s==0);
    
    % P = Un * Sn * Un', Un=U*diag(teken)
    [nr, ~] = size(U.'.*V');
    X = [zeros(nr,1) U.'.*V'];
    teken = sum(X')'; %#ok<UDIM>
    U = U .* teken(:,ones(1,n)).';
    teken(if0) = zeros(length(if0), 1);
    ifp = find(teken>0);
    ifn = find(teken<0);
    
    if ~isempty(ifp)
        Up = U(:, ifp);
        sp = sqrt(s(ifp));
        xo_pp = sp(:, ones(1,n)).'.*Up;
    end
    
    if ~isempty(ifn)
        Un = U(:, ifn);
        sn = sqrt(s(ifn));
        xo_pn = sn(:, ones(1,n)).'.*Un;
    end
end
