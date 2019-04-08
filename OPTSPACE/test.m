%[X, S, Y, dist] = OptSpace1(Dist.^2, 5, 1000, 1e-4);
%D_a = X*S*Y';
D_a = IncOPTSPACE_Cow(Dist.^2, 5);
M_copy = D_a;
D_a = sqrt(D_a);
D_a = D_a - diag(diag(D_a));
error = norm(D_a - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
% Gram_a = GramMatrix(D_a);
% Gram = GramMatrix(Dist_Truth);
% error_G = norm(Gram_a - Gram, 'fro')/norm(Gram, 'fro');
coor = MDS(D_a, 3);
ViewMesh(coor, trg)
M_opt_1 = D_a;
error_opt_1 = error;
coor_opt_1 = coor;
save 1sampledCow coor_opt_1 error_opt_1 M_opt_1 M_copy
%ViewMesh(coor, trg)
%save 3sampledCow coor error D_a

% coor = smacof(Dist, 3);
% D_a = DistMatrix(coor);
% D_a = D_a - diag(diag(D_a));
% error = norm(D_a - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
% Gram_a = GramMatrix(D_a);
% Gram = GramMatrix(Dist_Truth);
% error_G = norm(Gram_a - Gram, 'fro')/norm(Gram, 'fro');
% ViewMesh(coor, trg)