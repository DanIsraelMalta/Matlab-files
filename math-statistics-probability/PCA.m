%{
% PCA returns a given matrix principle components, the principle components Z-score,
% the covariance matrix eigenvalues and T-squared statistics for each covariance
% matrix eigenvalue.
%
% xi_mat - data matrix
% xo_pca - data matrix principle components
% xo_z_score - principle components Z-score
% xo_cov_eig - data matrix covariance eigenvalues
% xo_t_square - covariance T-squared statistics
%
% Dan I. Malta 2016
%}
function [xo_pca, xo_z_score, xo_cov_eig, xo_t_square] = PCA(xi_mat)
    % housekeeping
    [m,n] = size(xi_mat);
    r = min(m - 1, n);     % max possible rank of xi_mat
    
    % xi_mat center and scale
    avg = mean(xi_mat);
    centerx = (xi_mat - avg(ones(m, 1), :));
    
    % SVD
    [U, xo_cov_eig, xo_pca] = svd(centerx ./ sqrt(m - 1), 0);
    
    % Z-score
    xo_z_score = centerx * xo_pca;
    
    % covariance eigenvalues
    xo_cov_eig = diag(xo_cov_eig) .^ 2;
    
    % if rank is smaller then number of columns - zero pad it
    if r < n
        xo_cov_eig = [xo_cov_eig(1 : r); zeros(n - r, 1)];
        xo_z_score(:, r + 1 : end) = 0;
    end
    
    % T-squared
    tmp = sqrt(diag(1 ./ xo_cov_eig(1 : r))) * xo_z_score(:, 1 : r)';
    xo_t_square = sum(tmp .* tmp)';
end
