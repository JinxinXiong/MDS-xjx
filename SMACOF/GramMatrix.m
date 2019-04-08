function Gram = GramMatrix(D)
%DΪһ�������
D2 = D.^2;  % D2������ƽ������
n = size(D, 1);  % n������������ĸ���
H = eye(n) - 1 / n * ones(n); % H: centering matrix with size n*n
Gram = -1/2 * H * D2 * H; % B��the gram matrix