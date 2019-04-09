%% ²ÎÊý
%               MaxIter     tol     tau     delta
% SMACOF        1000        1E-4    5*n     1.2*(n^2/m)
% SVT           1000        1E-4
% OPTSPACE      1000        1E-4
%%
clear
load 1KSphereUniformSampled
Num = 5;
%% 20%sampled
error1 = 0;
error2 = 0;
error3 = 0;

time1 = 0;
time2 = 0;
time3 = 0;
run GenerateDistance20Sphere

for i = 1: Num
    t1 = clock;
    coor1 = smacof(Dist, 3);
    t2 = clock;
    time1 = time1 + etime(t2, t1);
    D1 = DistMatrix(coor1);
    D1 = D1 - diag(diag(D1));
    error = norm(D1 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error1 = error1 + error;
    
    t1 = clock;
    D2 = SVT(Dist.^2);
    D2 = sqrt(D2);
    D2 = D2 - diag(diag(D2));
    coor2 = MDS(D2, 3);
    t2 = clock;
    time2 = time2 + etime(t2, t1);
    error = norm(D2 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error2 = error2 + error;
    
    t1 = clock;
    D3 = OPTSPACE(Dist.^2, 4);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    coor3 = MDS(D3, 3);
    t2 = clock;
    time3 = time3 + etime(t2, t1);
    error = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error3 = error3 + error;
    
    if(i == 1)
        save SphereCoor20 coor1 coor2 coor3 D1 D2 D3
    end
end

error1 = error1/Num;
error2 = error2/Num;
error3 = error3/Num;

time1 = time1/Num;
time2 = time2/Num;
time3 = time3/Num;

save SphereResult20 error1 error2 error3 time1 time2 time3
%% 10%sampled
error1 = 0;
error2 = 0;
error3 = 0;

time1 = 0;
time2 = 0;
time3 = 0;
run GenerateDistance10Sphere
for i = 1: Num
    t1 = clock;
    coor1 = smacof(Dist, 3);
    t2 = clock;
    time1 = time1 + etime(t2, t1);
    D1 = DistMatrix(coor1);
    D1 = D1 - diag(diag(D1));
    error = norm(D1 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error1 = error1 + error;
    
    t1 = clock;
    D2 = SVT(Dist.^2);
    D2 = sqrt(D2);
    D2 = D2 - diag(diag(D2));
    coor2 = MDS(D2, 3);
    t2 = clock;
    time2 = time2 + etime(t2, t1);
    error = norm(D2 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error2 = error2 + error;
    
    t1 = clock;
    D3 = OPTSPACE(Dist.^2, 4);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    coor3 = MDS(D3, 3);
    t2 = clock;
    time3 = time3 + etime(t2, t1);
    error = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error3 = error3 + error;
    
    if(i == 1)
        save SphereCoor10 coor1 coor2 coor3 D1 D2 D3
    end
end

error1 = error1/Num;
error2 = error2/Num;
error3 = error3/Num;

time1 = time1/Num;
time2 = time2/Num;
time3 = time3/Num;

save SphereResult10 error1 error2 error3 time1 time2 time3
%% 5%sampled
error1 = 0;
error2 = 0;
error3 = 0;

time1 = 0;
time2 = 0;
time3 = 0;
run GenerateDistance5Sphere
for i = 1: Num
    t1 = clock;
    coor1 = smacof(Dist, 3);
    t2 = clock;
    time1 = time1 + etime(t2, t1);
    D1 = DistMatrix(coor1);
    D1 = D1 - diag(diag(D1));
    error = norm(D1 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error1 = error1 + error;
    
    t1 = clock;
    D2 = SVT(Dist.^2);
    D2 = sqrt(D2);
    D2 = D2 - diag(diag(D2));
    coor2 = MDS(D2, 3);
    t2 = clock;
    time2 = time2 + etime(t2, t1);
    error = norm(D2 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error2 = error2 + error;
    
    t1 = clock;
    D3 = OPTSPACE(Dist.^2, 4);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    coor3 = MDS(D3, 3);
    t2 = clock;
    time3 = time3 + etime(t2, t1);
    error = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error3 = error3 + error;
    
    if(i == 1)
        save SphereCoor5 coor1 coor2 coor3 D1 D2 D3
    end
end

error1 = error1/Num;
error2 = error2/Num;
error3 = error3/Num;

time1 = time1/Num;
time2 = time2/Num;
time3 = time3/Num;

save SphereResult5 error1 error2 error3 time1 time2 time3
%% 3%sampled
error1 = 0;
error2 = 0;
error3 = 0;

time1 = 0;
time2 = 0;
time3 = 0;
run GenerateDistance3Sphere
for i = 1: Num
    t1 = clock;
    coor1 = smacof(Dist, 3);
    t2 = clock;
    time1 = time1 + etime(t2, t1);
    D1 = DistMatrix(coor1);
    D1 = D1 - diag(diag(D1));
    error = norm(D1 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error1 = error1 + error;
    
    t1 = clock;
    D2 = SVT(Dist.^2);
    D2 = sqrt(D2);
    D2 = D2 - diag(diag(D2));
    coor2 = MDS(D2, 3);
    t2 = clock;
    time2 = time2 + etime(t2, t1);
    error = norm(D2 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error2 = error2 + error;
    
    t1 = clock;
    D3 = OPTSPACE(Dist.^2, 4);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    coor3 = MDS(D3, 3);
    t2 = clock;
    time3 = time3 + etime(t2, t1);
    error = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error3 = error3 + error;
    
    if(i == 1)
        save SphereCoor3 coor1 coor2 coor3 D1 D2 D3
    end
end

error1 = error1/Num;
error2 = error2/Num;
error3 = error3/Num;

time1 = time1/Num;
time2 = time2/Num;
time3 = time3/Num;

save SphereResult3 error1 error2 error3 time1 time2 time3
%% 2%sampled
error1 = 0;
error2 = 0;
error3 = 0;

time1 = 0;
time2 = 0;
time3 = 0;
run GenerateDistance2Sphere
for i = 1: Num
    t1 = clock;
    coor1 = smacof(Dist, 3);
    t2 = clock;
    time1 = time1 + etime(t2, t1);
    D1 = DistMatrix(coor1);
    D1 = D1 - diag(diag(D1));
    error = norm(D1 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error1 = error1 + error;
    
    t1 = clock;
    D2 = SVT(Dist.^2);
    D2 = sqrt(D2);
    D2 = D2 - diag(diag(D2));
    coor2 = MDS(D2, 3);
    t2 = clock;
    time2 = time2 + etime(t2, t1);
    error = norm(D2 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error2 = error2 + error;
    
    t1 = clock;
    D3 = OPTSPACE(Dist.^2, 4);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    coor3 = MDS(D3, 3);
    t2 = clock;
    time3 = time3 + etime(t2, t1);
    error = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error3 = error3 + error;
    
    if(i == 1)
        save SphereCoor2 coor1 coor2 coor3 D1 D2 D3
    end
end

error1 = error1/Num;
error2 = error2/Num;
error3 = error3/Num;

time1 = time1/Num;
time2 = time2/Num;
time3 = time3/Num;

save SphereResult2 error1 error2 error3 time1 time2 time3
