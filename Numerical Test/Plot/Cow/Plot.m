load CowError
plot(error1(1: 500), 'LineWidth', 1)
hold on 
plot(error2(1: 500), 'LineWidth', 1)
plot(error3(1: 500), 'LineWidth', 1)

%title("�����㷨��Cow���ݼ���Error����������仯����")

xlabel("��������")
ylabel("Error")
legend("SMACOF", "SVT", "OPTSPACE")