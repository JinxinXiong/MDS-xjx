load PlotError

plot(error1, 'LineWidth', 1)
hold on 
plot(error2, 'LineWidth', 1)
plot(error3, 'LineWidth', 1)

%title("�����㷨��S^2���ݼ�Error����������仯����")

xlabel("��������")
ylabel("Error")
legend("SMACOF", "SVT", "OPTSPACE")
