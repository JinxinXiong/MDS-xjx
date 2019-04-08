clear
load 1KSphereUniformSampled
run GenerateDistance5Sphere

MaxIter = 500;
epsilon = 1e-8;
tol = 1e-8;
error1 = zeros(1, MaxIter);
error2 = zeros(1, MaxIter);
error3 = zeros(1, MaxIter);
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
    B_Z(isnan(B_Z)) = 0;   % �������ڵ����ת����0
    B_Z = B_Z - diag(diag(B_Z)); % ���Խ�Ԫ�ر��0
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
r = 4;

[m, n] = size(M);
E = logical(M);
M_E = M;

num = sum(sum(M_E)); % �����е���degree

avg_row = 2 * num/ m;
avg_col = 2 * num/ n;

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

%S_0 = S(1:r, 1:r)/eps;  %�����в���Ҫ

%MaxIter = 1000;
%tol = 1e-4;
tau = 1e-3; % ������ֵ

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
%%
function S = optimize_F_S(X, Y, M_E, E)
% X, Y: the input value of variable in the function F
% M_E�� the observed matrix
% E: the mask matrix
% S��the matrix S that minimize the function F given value X, Y
% ���Ż���A*S_ij = b����ʽ���

[~, k] = size(X);
b = X'* M_E * Y;
b = b(:); %����������

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
%������ʧ����
F = sum( sum( ( (X*S*Y' - M_E).*E ).^2 ) )/2 ;
end