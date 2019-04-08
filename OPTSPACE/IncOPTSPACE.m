function M_approx = IncOPTSPACE(M,r_max)
%输入：
%   M: 待恢复的矩阵the observed matrix
%   r_max：即对秩循环的最大值，等同于希望恢复的结果的秩
%输出：
%   M_approx: 恢复出的矩阵
% 
% Cow 数据集中建议参数，手动调参结果，供参考
% sampling rate:    2%      3%      5%      10%     20%      
% trimming:         1       1       0.8     0.7     0.8
% tau:              1e-2    1e-2    1e-3    1e-3    1e-3
%Trim
[m, n] = size(M);
E = logical(M);
M_E = M;
num = sum(sum(M_E)); % 矩阵中的总degree

% avg_row = 2 * num/ m;
% avg_col = 2 * num/ n;

%cow中使用1倍的average
avg_row = 0.7 * num/ m;
avg_col = 0.7 * num/ n;

M_row = sum(M_E, 2);
M_col = sum(M_E, 1);
index_row = M_row > avg_row;
index_col = M_col > avg_col;
M_trim = M;
M_trim(index_row,:) = 0;
M_trim(:, index_col) = 0;

M_approx = zeros(size(M));
%%
for rho = 0: r_max
    %compute th rank-1 projection of Tr(M) - M_approx, as the initial value
    [X_0, S, Y_0] = svds(M_trim - M_approx, 1);
    X_0 = sqrt(m) * X_0;
    Y_0 = sqrt(n) * Y_0;
    
    if(rho ~= 0)
        X = [X, X_0];
        Y = [Y, Y_0];
    else
        X = X_0;
        Y = Y_0;
    end
    
    MaxIter = 500;
    tol = 1e-4;
    tau = 1e-1; % value to initialize the step

    for k = 1: MaxIter
        S = optimize_F_S(X, Y, M_E, E);
    
        grad_X = (E.* (X * S * Y' - M_E)) * Y * S';
        grad_Y = (E.* (X * S * Y' - M_E)') * X * S;
    
        t = tau;
    
        for i = 1: 50
            if(cost_F(X - t*grad_X, Y - t*grad_Y,S,M_E, E) - cost_F(X, Y,S,M_E, E) < ...
                t * 0.5 * (norm(grad_X, 'fro')^2 + norm(grad_Y, 'fro')^2))
                break;
            end
            t = t/2;
        end
        
        X_old = X;
        Y_old = Y;
        
        X = X - t*grad_X;
        Y = Y - t*grad_Y;
        eps = norm(E.*(M_E - X*S*Y'), 'fro') / norm(M_E, 'fro');
        %eps = norm(E.*(M_E - X*S*Y'), 'fro') / norm(M_E, 'fro');
        fprintf('rho = %d, %dth iteration, cost_F = %f,eps = %f, differ = %f\n', rho,k,cost_F(X_old, Y_old,S,M_E, E),eps, abs(cost_F(X,Y,S,M_E, E) - cost_F(X_old, Y_old,S,M_E, E)));
        if (1e-6 * cost_F(X_old, Y_old,S,M_E, E) - abs(cost_F(X,Y,S,M_E, E) - cost_F(X_old, Y_old,S,M_E, E)) > 1e-5 || eps < 1e-5)
            break;
        end
        %M_approx = X * S * Y';
        
        %fprintf('rho = %d, %dth iteration, eps: %f\n', rho, k, eps);
        
%         if(eps < tol)
%             fprintf('rho = %d convergent\n',rho);
%             break;
%         end
    end
    
    M_approx = X * S * Y';
    eps = norm(E.*(M_E - X*S*Y'), 'fro') / norm(M_E, 'fro');
    
    fprintf('rho = %d, %f \n', rho,eps);
    
    if(eps < tol)
        fprintf('rho = %d convergent\n',rho);
        break;
    else
        fprintf('Not convergent\n');
    end
end
end
%%
function S = optimize_F_S(X, Y, M_E, E)
% X, Y: the input value of variable in the function F
% M_E： the observed matrix
%E: the mask matrix
% return S, the matrix S that minimize the function F given value X, Y
% 想着化成A*S_ij = b的形式求解

[~, k] = size(X);
b = X'* M_E * Y;
b = b(:); %拉成列向量

A = zeros(k*k, k*k);
for i = 1: k
    for j = 1:k
        index = (j-1) * k + i;
        coef = X' * (X(:,i) * Y(:, j)' .* E) * Y;
        A(:, index) = coef(:);
    end
end
S = A\b;
S = reshape(S, k, k);
end
%%
function F = cost_F(X, Y, S, M_E, E)
F = sum( sum( ( (X*S*Y' - M_E).*E ).^2 ) )/2 ;
end









