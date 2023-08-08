clear;
close all;
clc;

%常数预设
a = 5.43 * 10 ^ (-10);
h = 6.63 * 10 ^ (-34);
m0 = 9.1 * 10 ^ (-31);
c = (h / 2 / pi) ^ 2 / 2 / m0;

syms x;
V = 1 * 10 ^ (-19) * cos(2 * pi / a * x);
n_size = 5; % n_size需要是奇数才能用特征根法
T_start = -a / 4;
T_end = a / 4;
T = a;
acc = 1 / 500;
p = 1 / 5;

eigs = my_eig(V, n_size, T_start, T_end, T, acc);
e = appro(V, n_size, T_start, T_end, T, acc, p);

%作图
figure;
Vx = subs(V, x, (-1.5 : 0.01 : 1.5) * a);
plot(-1.5 : 0.01 : 1.5, Vx);
xlabel('x/a');
ylabel('V(x)/J');
ylim([0, 1.5 * 10 ^ (-19)]);
title('势能曲线');

figure;
hold on;
for i = 1 : n_size
    plot(-1 : acc : 1, eigs(i, :));
end
xlabel('k(\pi/a)');
ylabel('E/J');
title('简约布里渊区图景下能带（特征根求解）');

figure;
hold on;
for i = 1 : n_size
    plot(-1 : acc : 1, e(i, :));
end
xlabel('k(\pi/a)');
ylabel('E/J');
title('简约布里渊区图景下能带（近自由电子近似求解）');

err = abs(e - eigs);
err_rela = err ./ eigs;

figure;
hold on;
for i = 1 : n_size
    plot(-1 : acc : 1, err(i, :), 'DisplayName', '第' + string(i) + '能带');
end
xlabel('k(\pi/a)');
ylabel('E/J');
title('绝对误差');
legend();

figure;
hold on;
for i = 1 : n_size
    plot(-1 : acc : 1, err_rela(i, :), 'DisplayName', '第' + string(i) + '能带');
end
xlabel('k(\pi/a)');
title('相对误差');
legend();

gap_eig = zeros(n_size - 1, 1);
disp('特征根法带隙:');
for i = 1 : n_size - 1
    gap_eig(i) = min(abs(eigs(i, :) - eigs(i + 1, :)));
    fprintf('%d: %d\n', i, gap_eig(i));
end

gap_appro = zeros(n_size - 1, 1);
disp('近自由电子近似法带隙:');
for i = 1 : n_size - 1
    gap_appro(i) = min(abs(e(i, :) - e(i + 1, :)));
    fprintf('%d: %d\n', i, gap_appro(i));
end
