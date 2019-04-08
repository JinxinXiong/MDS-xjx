load CowError
plot(error1(1: 500), 'LineWidth', 1)
hold on 
plot(error2(1: 500), 'LineWidth', 1)
plot(error3(1: 500), 'LineWidth', 1)

%title("三种算法在Cow数据集上Error随迭代次数变化曲线")

xlabel("迭代次数")
ylabel("Error")
legend("SMACOF", "SVT", "OPTSPACE")