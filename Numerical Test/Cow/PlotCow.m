load cowPointCloud
run GenerateDistance2Sphere
load CowCoor2
figure('color',[1 1 1]);
subplot(4, 4, 1)
ViewMesh(coor1, trg)

subplot(4, 4, 5)
ViewMesh(coor2, trg)

subplot(4, 4, 9)
ViewMesh(coor3, trg)

subplot(4, 4, 13)
ViewMesh(coor4, trg)

run GenerateDistance3Sphere
load CowCoor3
subplot(4, 4, 2)
ViewMesh(coor1, trg)

subplot(4, 4, 6)
ViewMesh(coor2, trg)

subplot(4, 4, 10)
ViewMesh(coor3, trg)

subplot(4, 4, 14)
ViewMesh(coor4, trg)

run GenerateDistance5Sphere
load CowCoor5
subplot(4, 4, 3)
ViewMesh(coor1, trg)

subplot(4, 4, 7)
ViewMesh(coor2, trg)

subplot(4, 4, 11)
ViewMesh(coor3, trg)

subplot(4, 4, 15)
ViewMesh(coor4, trg)

run GenerateDistance10Sphere
load CowCoor10
subplot(4, 4, 4)
ViewMesh(coor1, trg)

subplot(4, 4, 8)
ViewMesh(coor2, trg)

subplot(4, 4, 12)
ViewMesh(coor3, trg)

subplot(4, 4, 16)
ViewMesh(coor4, trg)