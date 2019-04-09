function [coor] = smacof(D, d)
% coor: return the approximate coordinate
% D: the input distance matrix, not need to be full 
% d: the dimension of the space intended to embed in

epsilon = 1e-4; 
MaxIter = 1000; 
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
end

if(k == MaxIter + 1)
    fprintf('not convergent');
else
    fprintf('convergent');
end
coor = X;
