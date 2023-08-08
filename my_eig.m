function eigs = my_eig(V, n_size, T_start, T_end, T, acc)
% my_eig:特征根法
% V:势能函数
% n_size:需要求解的能带数目
% T_start:计算傅里叶级数积分的起始坐标
% T_end:计算傅里叶级数积分的终止坐标
% T:势能函数的周期
% acc:希望获得的计算结果中k的精度（以pi/T为单位）

%常数预设
h = 6.63 * 10 ^ (-34);
m0 = 9.1 * 10 ^ (-31);
c = (h / 2 / pi) ^ 2 / 2 / m0;

%求解Vn
syms x n;
% 含n的表达式，int是在做积分
Vn = 1 / T * int(V * exp(-1i * n * 2 * pi * x / T), x, T_start, T_end);
%仅求解[-n_size,n_size]之间的傅里叶级数
Vnumber = zeros(2 * n_size + 1, 1);
for k = -n_size : n_size
    Vnumber(k + n_size + 1) = double(subs(Vn, n, k)); 
end

%滤除数值计算带来的极小值（解析求解为0值）
maxVn=max(abs(Vnumber));
Vnumber(abs(Vnumber)<(maxVn*10^(-10)))=0;

%在-pi/T到pi/T之间求解k的本征能量
k_sery = -pi / T : acc * pi / T : pi / T;
l = length(k_sery);
%等效为求解特征值
eigs = zeros(n_size,l);
for i = 1 : l
    k = k_sery(i);
    M = zeros(n_size, n_size);
    % 为M赋值
    for j = 1 : n_size
        M(j, j) = c * (k + 2 * pi / T * (j - (n_size + 1) / 2)) ^ 2 ...
        + Vnumber(n_size + 1);
    end
    for j = 1 : n_size - 1
        for k = 1 : n_size - j
            M(k, k + j) = Vnumber(n_size + 1 - j);
            M(k + j, k) = Vnumber(n_size + 1 + j);
        end
    end
    % eig可以求解特征值
    eigs(:, i) = eig(M);
end

end