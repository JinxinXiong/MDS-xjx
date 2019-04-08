function Gram = GramMatrix(D)
%D为一距离矩阵
D2 = D.^2;  % D2：距离平方矩阵
n = size(D, 1);  % n：列数，即点的个数
H = eye(n) - 1 / n * ones(n); % H: centering matrix with size n*n
Gram = -1/2 * H * D2 * H; % B：the gram matrix