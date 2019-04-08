load PlotError

plot(error1, 'LineWidth', 1)
hold on 
plot(error2, 'LineWidth', 1)
plot(error3, 'LineWidth', 1)

%title("三种算法在S^2数据集Error随迭代次数变化曲线")

xlabel("迭代次数")
ylabel("Error")
legend("SMACOF", "SVT", "OPTSPACE")
