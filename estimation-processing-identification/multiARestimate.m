%{
% multiARestimate perform multivariate (vector) adaptive AR process
% estimation based upon a multidimensional kalman filter algorithm.
%
% the standard multivariate AR process is defines as:
% y[n] - A1[n] * y[n-1] - ... - Ap[n] * y[n-p] = e[n]
% where y[n] dimension is s.
% 
% y - signal (process)
% p - maximal model order
% uc - update coefficient (0 < uc << 1)
% mode - process noise update method [0, 4]
% x - process state vector (dimension is s * s * p)
% e - predicition error (dimension is s)
% Kalman - structure holding latest kalman estimator (x & Q2 are derived from it)
% Q2 - estimated measurement noise covariance matrix (dimension sxs)
%
% Dan I. Malta 2013
%}
function [x, e, Kalman, Q2] = multiARestimate(y, p, UC, mode)
    % housekeeping
    [M, LEN] = size(y');		%number of channels, total signal length
    L = M * M * p;
    if LEN<(p+1),
        disp('Not enough observed data supplied for given model order.\n');
        return
    end

    % initialization
    ye = zeros(size(y));	% output predictor
    x=zeros(L,LEN);
    Q2=zeros(M, M, LEN);
    upd = eye(L) / L * UC;		%diagonal matrix containing UC
    
    % kalman filter initialization
    Kalman = struct('F', eye(L),...
                                          'H', zeros(M, L),...
                                          'G', zeros(L, M),...
                                          'x', zeros(L, 1),...
                                          'Kp', eye(L),...  % K predicted or a-priori given = K(n+1, n)
                                          'Q1', eye(L) * UC,...
                                          'Q2', eye(M),...
                                          'ye', zeros(M, 1));

    % 3/4 modes
    if mode == 3
        Block = kron(eye(M), ones(M * p));
    elseif mode == 4
        index = [];
        Block1 = [];
        Block0 = [];
        for i = 1 : M
                index = [index ((i - 1) * M * p + (i : M : i * M * p))];
                mone = eye(M);
                mone(i, i) = 0;
                mzero = eye(M) - mone;
                Block1 = Blkdiag(Block1, kron(eye(p), mone));
                Block0 = Blkdiag(Block0, kron(eye(p), mzero));
        end
    end

    % multivariate kalman estimation
    for n = 2 : LEN
        % past observations
        if n <= p
            Yr = [y(n - 1 : -1 : 1, :)' zeros(M, p - n + 1)];
            Yr = Yr(:)';
        else
            Yr = y(n - 1 : -1 : n - p, :)';
            Yr = Yr(:)';
        end
        
        % measurement matrix
        Kalman.H = kron(eye(M), Yr);
        
        % prediction error
        ye(n, :) = (Kalman.H * Kalman.x)';
        err=y(n, :) - ye(n, :);
        
        % NaN's are trouble in this kind of schema
        if ~any(isnan(err(:))),
                %update of Q2 using the prediction error of the previous step
                Kalman.Q2 = (1 - UC) * Kalman.Q2 + UC * err' * err;
                
                % a-posteriori state error covariance matrix
                KpH = Kalman.Kp * Kalman.H';
                HKp = Kalman.H * Kalman.Kp;
                K = Kalman.Kp - Kalman.G * HKp; 
                
                %Kalman gain
                Kalman.G = KpH * inv(Kalman.H * KpH + Kalman.Q2);
                
                % Q1 update using the predicted state error cov matrix (in mode 0 Q1 is not updated)
                if mode == 1
                    Kalman.Q1 = diag(diag(K)) .* UC;
                elseif mode == 2
                        Kalman.Q1 = upd * trace(K);
                elseif mode == 3
                        Kalman.Q1 = diag(sum((Block * diag(diag(K)))')) / (p * M) * UC; %#ok<UDIM>
                elseif mode == 4
                        avg = trace(K(index, index)) / (p * M) * UC;
                        Kalman.Q1 = Block1 * UC + Block0 * avg;
                end
                
                % a-priori state error covariance matrix for the next time step
                Kalman.Kp = K + Kalman.Q1;
                
                % current estimation of state x
                Kalman.x = Kalman.x + Kalman.G * (err)';
        end
        
        % output extraction
        x(:,n) = Kalman.x;
        Q2(:,:,n) = Kalman.Q2;
    end
    
    e = y - ye;
    x = x';
end
