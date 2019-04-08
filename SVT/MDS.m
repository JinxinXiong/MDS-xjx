function [coor] = MDS(D, p)
%description: 函数MDS用MDS算法，通过距离矩阵来恢复坐标
%传入参数：
%   D：距离矩阵
%   p：需要返回的坐标空间的维数
%输出参数：
%   coor：坐标点坐标矩阵
%--------------------------------------------------------------------------

%--------------------------------SETUP-------------------------------------
D2 = D.^2;  % D2：距离平方矩阵
n = size(D, 1);  % n：列数，即点的个数
H = eye(n) - 1 / n * ones(n); % H: centering matrix with size n*n
B = -1/2 * H * D2 * H; % B：the gram matrix

[eigvec,eigval] = eig(B);
% eigval：每一列为B的特征向量
% eigval：为一对角矩阵，每一个元素为一个特征值
[eigval_sort,index] = sort(diag(eigval),'descend');
eigval_sort = diag(eigval_sort(index));
eigvec_sort = eigvec(:,index);

Eigvec = eigvec_sort(:, 1:p);
Eigval = eigval_sort(1:p, 1:p);

coor = Eigvec * sqrt(Eigval);
end

