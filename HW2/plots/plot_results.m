figure;
T = readtable("out_Gaussian.csv");
labels = T.Properties.VariableNames;
A = table2array(T);
subplot(2,1,1);
plot(A(:,1), [A(:,2),A(:,4), A(:,6)]);
title("Execution time");
xlabel(labels(1));
ylabel("Time");
legend(labels([2,4,6]));
subplot(2,1,2);
plot(A(:,1), [A(:,3),A(:,5), A(:,7)]);
title("Error");
xlabel(labels(1));
ylabel("Ax-b");
legend(labels([3,5,7]));

savefig("Gaussian_plot.fig");


figure;
T = readtable("out_Cholesky.csv");
labels = T.Properties.VariableNames;
A = table2array(T);
subplot(2,1,1);
plot(A(:,1), [A(:,2),A(:,4)]);
title("Execution time");
xlabel(labels(1));
ylabel("Time");
legend(labels([2,4]));
subplot(2,1,2);
plot(A(:,1), [A(:,3),A(:,5)]);
title("Error");
xlabel(labels(1));
ylabel("Ax-b");
legend(labels([3,5]));



savefig("Cholesky_plot.fig");