clc; clear;
load('observation.mat', 'pts_o');
load('pts.mat', 'pts_markers');
load('dist.mat', 'dist');
load('gt.mat', 'pts_marks_gt')

lambda = 0.0168318035;

%% ---------------- Apply Newton Optimization ----------------- 
p_opt = optimize(pts_o, pts_markers, dist, lambda, 150);
rmse = calculate_rmse(p_opt, pts_marks_gt);
fprintf('Lambda = %.10f, RMSE = %.10f\n', lambda, rmse);

%% ----------------- Calculate optimal lambda ------------------ 
% Uncomment this section to find the optimal lambda
% lambda_values = logspace(-5, 1, 200);
% rmse_values = zeros(size(lambda_values));
%     
% for j = 1:length(lambda_values)
%     p_opt = optimize(pts_o, pts_markers, dist, lambda_values(j), 150);
%     rmse = calculate_rmse(p_opt, pts_marks_gt);
%     rmse_values(j) = rmse;
% end
%     
% semilogx(lambda_values, rmse_values, '.-');
% title('RMSE vs \lambda');
% xlabel('\lambda');
% ylabel('RMSE');
% 
% [~, index] = min(rmse_values);
% lambda = lambda_values(index);
% fprintf('Minimum lambda, lambda_op = %.10f\n', lambda);

%% ----------------- Plot RMSE vs Lambda ------------------
% Uncomment this section to compare different initialization schemes
% k = 1:150;
% rmse = zeros(length(k),1);
% for k = k
%     p_opt = optimize(pts_o, pts_markers, dist, lambda, k);
%     rmse(k) = calculate_rmse(p_opt, pts_marks_gt);
% end 
% figure;
% plot((1:150), rmse);
% title('Number of iterations vs RMSE');
% xlabel('Number of iterations'); ylabel('RMSE');
% xticks(0:25:150); yticks(0:21);

%% Newton's Optimization
function p_opt = optimize(pts_o, pts_markers, dist, lambda, max_iter)
    alpha = 0.1;
        
    % Number of markers on tunnel surface
    N = size(pts_markers, 2);
    % Number of LiDAR markers on ground    
    K = size(pts_o, 1);
        
    % Initialize matrix to store estimated marker positions
    p_opt = zeros(N, 3);
    for i = 1:N
        % Noisy position of the i-th marker seen by K LiDARs        
        p_hat = squeeze(pts_markers(:, i, :));
        
        % Distance between i-th marker and K LiDARs         
        d = dist(i, :)';
                
        % Initialize marker position as average of the noisy coordinates
        p = mean(p_hat);
        
        % Initialize marker as single coordinate from the noise measurements                
%         p = p_hat(1,:);
        
        % Initialize marker position as random vector              
%         p = rand(1,3);
            
        % Keep track of iterations of Newton's method
        iter = 0;
                
        % While Newton's method is not converged apply Newton's method
        while iter < max_iter
            % Vector of residuals
            r = vecnorm(p - pts_o, 2, 2) - d;
            % Jacobian
            J = (p - pts_o) ./ vecnorm(p - pts_o, 2, 2);
                        
            % Gradient of loss function/error function
            g = J' * r + lambda * sum(2 * (p - p_hat))';
            % Hessian of loss function
            H = J' * J + 2 * lambda * K * eye(3);
                        
            % Descent direction for Newton’s method
            delta_p = -inv(H)*g;
            % Updated iteration of p
            p = p + alpha * delta_p';
            % Place optimized coordinates in matrix
            p_opt(i, :) = p;
            
            iter = iter + 1;
            
            if norm(delta_p) < 1e-6
	            break;
            end
        end
    end
end

function rmse = calculate_rmse(p_opt, p)
    rmse = sqrt(norm(p_opt - p, 'fro').^2/length(p));
end