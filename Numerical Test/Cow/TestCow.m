clear
load cowPointCloud
Num = 5;
error1 = 0;
error2 = 0;
error3 = 0;
error4 = 0;

time1 = 0;
time2 = 0;
time3 = 0;
time4 = 0;

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
    D3 = OPTSPACE20(Dist.^2, 5);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    coor3 = MDS(D3, 3);
    t2 = clock;
    time3 = time3 + etime(t2, t1);
    error = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error3 = error3 + error;
    
    t1 = clock;
    D4 = IncOPTSPACE20(Dist.^2, 5);
    D4 = sqrt(D4);
    D4 = D4 - diag(diag(D4));
    coor4 = MDS(D4, 3);
    t2 = clock;
    time4 = time4 + etime(t2, t1);
    error = norm(D4 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error4 = error4 + error;
    
    if(i == 1)
        save CowCoor20 coor1 coor2 coor3 coor4 D1 D2 D3 D4
    end
end

error1 = error1/Num;
error2 = error2/Num;
error3 = error3/Num;
error4 = error4/Num;

time1 = time1/Num;
time2 = time2/Num;
time3 = time3/Num;
time4 = time4/Num;

save CowResult20 error1 error2 error3 error4 time1 time2 time3 time4
%%
error1 = 0;
error2 = 0;
error3 = 0;
error4 = 0;

time1 = 0;
time2 = 0;
time3 = 0;
time4 = 0;
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
    D3 = OPTSPACE10(Dist.^2, 5);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    coor3 = MDS(D3, 3);
    t2 = clock;
    time3 = time3 + etime(t2, t1);
    error = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error3 = error3 + error;
    
    t1 = clock;
    D4 = IncOPTSPACE10(Dist.^2, 5);
    D4 = sqrt(D4);
    D4 = D4 - diag(diag(D4));
    coor4 = MDS(D4, 3);
    t2 = clock;
    time4 = time4 + etime(t2, t1);
    error = norm(D4 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error4 = error4 + error;
    
    if(i == 1)
        save CowCoor10 coor1 coor2 coor3 coor4 D1 D2 D3 D4
    end
end

error1 = error1/Num;
error2 = error2/Num;
error3 = error3/Num;
error4 = error4/Num;

time1 = time1/Num;
time2 = time2/Num;
time3 = time3/Num;
time4 = time4/Num;
save CowResult10 error1 error2 error3 error4 time1 time2 time3 time4
%%
error1 = 0;
error2 = 0;
error3 = 0;
error4 = 0;

time1 = 0;
time2 = 0;
time3 = 0;
time4 = 0;
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
    D3 = OPTSPACE5(Dist.^2, 5);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    coor3 = MDS(D3, 3);
    t2 = clock;
    time3 = time3 + etime(t2, t1);
    error = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error3 = error3 + error;
    
    t1 = clock;
    D4 = IncOPTSPACE5(Dist.^2, 5);
    D4 = sqrt(D4);
    D4 = D4 - diag(diag(D4));
    coor4 = MDS(D4, 3);
    t2 = clock;
    time4 = time4 + etime(t2, t1);
    error = norm(D4 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error4 = error4 + error;
    
    if(i == 1)
        save CowCoor5 coor1 coor2 coor3 coor4 D1 D2 D3 D4
    end
end

error1 = error1/Num;
error2 = error2/Num;
error3 = error3/Num;
error4 = error4/Num;

time1 = time1/Num;
time2 = time2/Num;
time3 = time3/Num;
time4 = time4/Num;


save CowResult5 error1 error2 error3 error4 time1 time2 time3 time4
%%
error1 = 0;
error2 = 0;
error3 = 0;
error4 = 0;

time1 = 0;
time2 = 0;
time3 = 0;
time4 = 0;
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
    D3 = OPTSPACE3(Dist.^2, 5);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    coor3 = MDS(D3, 3);
    t2 = clock;
    time3 = time3 + etime(t2, t1);
    error = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error3 = error3 + error;
    
    t1 = clock;
    D4 = IncOPTSPACE3(Dist.^2, 5);
    D4 = sqrt(D4);
    D4 = D4 - diag(diag(D4));
    coor4 = MDS(D4, 3);
    t2 = clock;
    time4 = time4 + etime(t2, t1);
    error = norm(D4 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error4 = error4 + error;
    
    if(i == 1)
        save CowCoor3 coor1 coor2 coor3 coor4 D1 D2 D3 D4
    end
end

error1 = error1/Num;
error2 = error2/Num;
error3 = error3/Num;
error4 = error4/Num;

time1 = time1/Num;
time2 = time2/Num;
time3 = time3/Num;
time4 = time4/Num;
save CowResult3 error1 error2 error3 error4 time1 time2 time3 time4
%%
error1 = 0;
error2 = 0;
error3 = 0;
error4 = 0;

time1 = 0;
time2 = 0;
time3 = 0;
time4 = 0;
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
    D3 = OPTSPACE2(Dist.^2, 5);
    D3 = sqrt(D3);
    D3 = D3 - diag(diag(D3));
    coor3 = MDS(D3, 3);
    t2 = clock;
    time3 = time3 + etime(t2, t1);
    error = norm(D3 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error3 = error3 + error;
    
    t1 = clock;
    D4 = IncOPTSPACE2(Dist.^2, 5);
    D4 = sqrt(D4);
    D4 = D4 - diag(diag(D4));
    coor4 = MDS(D4, 3);
    t2 = clock;
    time4 = time4 + etime(t2, t1);
    error = norm(D4 - Dist_Truth, 'fro')/norm(Dist_Truth, 'fro');
    error4 = error4 + error;
    
    if(i == 1)
        save CowCoor2 coor1 coor2 coor3 coor4 D1 D2 D3 D4
    end
end

error1 = error1/Num;
error2 = error2/Num;
error3 = error3/Num;
error4 = error4/Num;

time1 = time1/Num;
time2 = time2/Num;
time3 = time3/Num;
time4 = time4/Num;
save CowResult2 error1 error2 error3 error4 time1 time2 time3 time4
