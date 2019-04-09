load 1KSphereUniformSampled
run GenerateDistance2Sphere
load SphereCoor2
figure('color',[1 1 1]);
subplot(3, 4, 1)
ViewMesh(coor1, trg)

subplot(3, 4, 5)
ViewMesh(coor2, trg)

subplot(3, 4, 9)
ViewMesh(coor3, trg)

run GenerateDistance3Sphere
load SphereCoor3
subplot(3, 4, 2)
ViewMesh(coor1, trg)

subplot(3, 4, 6)
ViewMesh(coor2, trg)

subplot(3, 4, 10)
ViewMesh(coor3, trg)

run GenerateDistance5Sphere
load SphereCoor5
subplot(3, 4, 3)
ViewMesh(coor1, trg)

subplot(3, 4, 7)
ViewMesh(coor2, trg)

subplot(3, 4, 11)
ViewMesh(coor3, trg)

run GenerateDistance10Sphere
load SphereCoor10
subplot(3, 4, 4)
ViewMesh(coor1, trg)

subplot(3, 4, 8)
ViewMesh(coor2, trg)

subplot(3, 4, 12)
ViewMesh(coor3, trg)