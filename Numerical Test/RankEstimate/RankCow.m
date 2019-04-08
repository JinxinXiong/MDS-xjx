clear
%由于只需要看大致效果，将算法中的最大迭代次数改为500
load cowPointCloud
error30 = zeros(1, 5);
error50 = zeros(1, 5);
error70 = zeros(1, 5);
run GenerateDistance30Sphere
for r = 1: 5
    D3 = OPTSPACE(Dist.^2, r);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    error30(r) = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
end

run GenerateDistance50Sphere
for r = 1: 5
    D3 = OPTSPACE(Dist.^2, r);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    error50(r) = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
end

run GenerateDistance70Sphere
for r = 1: 5
    D3 = OPTSPACE(Dist.^2, r);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    error70(r) = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
end

save RankCow error30 error50 error70