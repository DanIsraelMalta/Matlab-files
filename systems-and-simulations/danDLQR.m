%{
% danDLQR is a variant to control toolbox function 'dlqr',
% as it performs the same operation but with the following options/modifications:
% 1) possibly to minimize either the riccaty equation (with a fudge parameter)
%    or a frequency weightd quadratic cost function.
% 2) enter the system in two different forms - either in the standard form (for control purposes),
%    or in the rational form (as retrieved from identification expermints).
%
% first form:
%    [Cfb, clP] = DOFLQ(P,H,lambda,Fy,Fu)
%    calculates the optimal output feedback transfer function such that the control law is:
%     u[k] = -Cfb(q) 
%     y[k] = P(q) * u[k] + H(q) * e[K] (e[K] is zero mean white noise)
%    minimises the frequency weighted cost function:
%     J = E {(Fy * y)^2 + lambda*(Fu * u)^2} (Fy, Fu - optional)
%    the output is:
%     Cfb(q) - optimal output
%     clP(q) - closed loop characteristic equation
%
% second form:
%    [S, R, clP] = DOFLQ(A,B,C,lambda)
%    calculates the optimal output feedback transfer function such that the control law is:
%     u[k] = -S(q)/R(q)
%     A(q)y[k] = B(q)u[k] + C(q)e[k]  (e[K] is zero mean white noise)
%    minimises the cost function
%     J = E {y^2 + lambda*u^2}
%    the output is:
%     cLp - closde loop characteristic equaiton
%     S - inifinite horizon solution associated with riccati
%     R - R matrix
% 
% Constraints: 
%  * P(q), H(q), Fy(q) and Fu(q) must be transfer functions.
%  * A(q), B(q) and C(q) must be polynomials described by row vectors with the same length.
%  * The first element of B(q) must be zero.
%  * The first element of A(q) cannot be zero.
%  * The polynomials cannot be zero.
%  * Lambda must be a non-negative scalar.
% 
%
% Dan I. Malta 2011
%}
function [out1, out2, out3] = danDLQR(par1, par2, par3, par4, par5)
  % first form (transfer function input)
  if isa(par1, 'tf')
    % controller frequency cost
    if nargin < 5
      Fu = 1;
    elseif isa(par5, 'tf')
      Fu = par5;
    end
    
    % output frequency cost
    if nargin < 4
      Fy = 1;
    elseif isa(par4, 'tf')
      Fy = par4;
    end
    
    % input
    lambda = par3;
    P = par1;
    H = par2;
   
    % minimal realization
    Pbar = minreal(Fy * inv(Fu) * P);
    Hbar = minreal(Fy * H);
    PHbar = minreal([ss(Pbar) ss(Hbar)]);
    
    % transfer data as vector
    [B, A] = tfdata(tf(PHbar(1,1)), 'v');
    [C, D] = tfdata(tf(PHbar(1,2)), 'v');
    
    % stabilize C
    rC = dsort(roots(C));
    while max(abs(rC)) > 1
      rC(1) = 1 / rC(1);
      rC = dsort(rC);
    end
    C = poly(rC);
   
    % call polynomial form
    [Slq, Rlq, PC] = danDLQR(A, B, C, lambda);
    out1 = minreal(tf(Slq, Rlq, -1) * Fy * inv(Fu));
    out2 = PC;
  else % first form (polynomial input)
    % housekeeping
    A = par1;
    B = par2;
    C = par3;
    lambda = par4;
    [mA, nA] = size(A); %#ok<ASGLU>
    [mB, nB] = size(B);
    [mC, nC] = size(C); %#ok<ASGLU>
   
    % A and C must be monic.
    if A(1) ~= 1
      C = C / A(1);
      B = B / A(1);
      A = A / A(1);
    end
    while C(1) == 0
      C = [C(2:nC) 0];
    end
    if C(1) ~= 1
      C = C / C(1);
    end
    
    % cost function
    rPP = lambda * conv(A, fliplr(A)) + conv(B, fliplr(B));
    rts = flipud(dsort(roots(rPP)));
    P = poly(rts(1 : nA-1));
    PC = conv(P, C);
    PP = conv(P, fliplr(P));
    r = sum(abs(rPP)) / sum(abs(PP));
    As = fliplr(A);
    rP = r * P;
    M = zeros(2*nA-2);
    for i = 1 : nA - 1
      M(i : i + nA - 1, i) = As';
      M(i : i + nA - 1, i + nA - 1) = rP';
    end
    BCs = conv(B(2 : nA), fliplr(C));
    XSs = M \ BCs';
    X(1 : nA) = [0 XSs(1 : nA - 1)'];
    Ss = [0 XSs(nA : 2 * nA - 2)'];
    S = real(fliplr(Ss));
   
    % solution
    B1 = B;
    while B1(1) == 0
      B1 = B1(2 : length(B1));
    end
    [Rs, res1] = deconv(conv(fliplr(P), X) + lambda * conv(A, Ss), B1);
    R = real(fliplr(Rs(length(Rs) - nA + 1 : length(Rs))));
    res2 = sum(abs(conv(A, R) + conv(B, S) - PC));
    if (sum(abs(res1)) > 1e-5) | (res2 > 1e-5) %#ok<OR2>
        disp('Warning: spectral factorization of residuals is too big. Soltion might be off.')
    end
   
    % output
    out1 = S;
    out2 = R;
    out3 = PC;
  end
end
