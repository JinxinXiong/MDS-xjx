function [coor] = MDS(D, p)
%description: ����MDS��MDS�㷨��ͨ������������ָ�����
%���������
%   D���������
%   p����Ҫ���ص�����ռ��ά��
%���������
%   coor��������������
%--------------------------------------------------------------------------

%--------------------------------SETUP-------------------------------------
D2 = D.^2;  % D2������ƽ������
n = size(D, 1);  % n������������ĸ���
H = eye(n) - 1 / n * ones(n); % H: centering matrix with size n*n
B = -1/2 * H * D2 * H; % B��the gram matrix

[eigvec,eigval] = eig(B);
% eigval��ÿһ��ΪB����������
% eigval��Ϊһ�ԽǾ���ÿһ��Ԫ��Ϊһ������ֵ
[eigval_sort,index] = sort(diag(eigval),'descend');
eigval_sort = diag(eigval_sort(index));
eigvec_sort = eigvec(:,index);

Eigvec = eigvec_sort(:, 1:p);
Eigval = eigval_sort(1:p, 1:p);

coor = Eigvec * sqrt(Eigval);
end

