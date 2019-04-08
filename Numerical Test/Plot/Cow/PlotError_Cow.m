clear
load cowPointCloud
run GenerateDistance5Sphere

MaxIter = 1000;
epsilon = 1e-8;
tol = 1e-8;
error1 = zeros(1, MaxIter);
error2 = zeros(1, MaxIter);
error3 = zeros(1, MaxIter);
error4 = zeros(1, 5 * MaxIter);
%% SMACOF
D = Dist;
d = 3;
%epsilon = 1e-4; 
%MaxIter = 1000; 
n = size(D,1); 
X = randn(n, d); 
Z = X;
W = logical(D) + eye(n); 
V = -W;
V = V + diag(sum(W, 2));
k = 0;
sigma = 1/2*(sum(sum(W.*(D - DistMatrix(X)).^2)));
sigma_old = zeros(size(sigma));
b = -W .* D; % used to calculate the matrix B_Z

while (k==0 ||(sigma_old - sigma > epsilon && k<= MaxIter))
    k = k + 1;
    sigma_old = sigma;
    
    B_Z = b./ DistMatrix(Z);
    B_Z(isnan(B_Z)) = 0;   % 将不存在的情况转换成0
    B_Z = B_Z - diag(diag(B_Z)); % 将对角元素变成0
    B_Z = B_Z + diag(- sum(B_Z, 2));
    
    X = pinv(V) * B_Z * Z;
    
    sigma = 1/2*(sum(sum(W.*(D - DistMatrix(X)).^2)));
    
    fprintf('%dth iteration: %e \n', k, sigma_old - sigma);
    Z = X;
    
    D1 = DistMatrix(X);
    D1 = D1 - diag(diag(D1));
    error1(k) = norm(D1 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
end
%% SVT
D = Dist.^2;

[n, ~] = size(D);
Omega = logical(D);
m = sum(sum(Omega));
%tau = 5 * n;
tau = 5 * n;   

delta = 1.2 * (n^2/m);
%delta = 10;

%MaxIter = 1000;
%tol = 1e-4;
l = 5;
k0 = ceil(tau/(delta*normest(D, 1e-2)));

Y = k0*delta*D;
r = 0;

for k = 1:MaxIter
    s = min([r + l, n]);
    
    flag = 0;
    while ~flag
        [U, S, V] = svds(Y, s);
        flag = (S(s,s) <= tau) || ( s == n );
        s = min(s + l, n);
    end
    sigma = diag(S); 
    r = sum(sigma > tau);
    
    X = U(:, 1:r) * diag(sigma(1:r)-tau)*V(:,1:r)';

    if(norm(Omega.*X - D, 'fro')/norm(D, 'fro') <= tol)
        fprintf('convergent');
        break;
    else
        fprintf('%dth iteration: %e\n',k,norm(Omega.*X - D, 'fro')/norm(D, 'fro'));
    end

    Y = Omega.*(Y + delta * (D - X));
    
    D2 = sqrt(X);
    D2 = D2 - diag(diag(D2));
    error2(k) = norm(D2 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
end
%% OPTSPACE
M = Dist.^2;
r = 5;

[m, n] = size(M);
E = logical(M);
M_E = M;

num = sum(sum(M_E)); % 矩阵中的总degree

avg_row = 0.8 * num/ m;
avg_col = 0.8 * num/ n;

M_row = sum(M_E, 2);
M_col = sum(M_E, 1);

index_row = M_row > avg_row;
index_col = M_col > avg_col;

M_trim = M;
M_trim(index_row,:) = 0;
M_trim(:, index_col) = 0;

[X, S, Y] = svds(M_trim, r);

X = sqrt(m) * X;
Y = sqrt(n) * Y;

%S_0 = S(1:r, 1:r)/eps;  %运算中不需要

%MaxIter = 1000;
%tol = 1e-4;
tau = 1e-2; % 步长初值

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
    
    X = X - t*grad_X;
    Y = Y - t*grad_Y;
    
    eps = norm(E.*(M_E - X*S*Y'), 'fro') / norm(M_E, 'fro');
    fprintf('%dth iteration: %e \n', k, eps);
    
    if (eps < tol)
        fprintf('convergent');
        break;
    end
    
    M_approx = X*S*Y';
    D3 = sqrt(M_approx);
    D3 = D3 - diag(diag(D3));
    error3(k) = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
end
%% IncOPTSPACE
M = Dist.^2;
r_max = 5;
[m, n] = size(M);
E = logical(M);
M_E = M;
num = sum(sum(M_E)); % 矩阵中的总degree

% avg_row = 2 * num/ m;
% avg_col = 2 * num/ n;

%cow中使用1倍的average
avg_row = 0.8 * num/ m;
avg_col = 0.8 * num/ n;

M_row = sum(M_E, 2);
M_col = sum(M_E, 1);
index_row = M_row > avg_row;
index_col = M_col > avg_col;
M_trim = M;
M_trim(index_row,:) = 0;
M_trim(:, index_col) = 0;

M_approx = zeros(size(M));
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
    
    MaxIter = 200;
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
        M_approx = X * S * Y';
        D4 = sqrt(M_approx);
        D4 = D4 - diag(diag(D4));
        error4(k) = norm(D4 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
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
    
    
end

%%
function S = optimize_F_S(X, Y, M_E, E)
% X, Y: the input value of variable in the function F
% M_E： the observed matrix
% E: the mask matrix
% S：the matrix S that minimize the function F given value X, Y
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

function F = cost_F(X, Y, S, M_E, E)
% 计算损失函数
F = sum( sum( ( (X*S*Y' - M_E).*E ).^2 ) )/2 ;
end