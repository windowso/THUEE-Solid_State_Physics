function e = appro(V, n_size, T_start, T_end, T, acc, p)
% appro:近似法
% V:势能函数
% n_size:需要求解的能带数目
% T_start:计算傅里叶级数积分的起始坐标
% T_end:计算傅里叶级数积分的终止坐标
% T:势能函数的周期
% acc:希望获得的计算结果中k的精度（以pi/T为单位）
% p:小于p认定为布里渊区“附近”

%常数预设
h = 6.63 * 10 ^ (-34);
m0 = 9.1 * 10 ^ (-31);
h_bar = h / (2 * pi);
n_size = n_size + 1; % 为了保持对称性，需要多计算一条能带

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

k = (-n_size : acc : n_size + 1) * pi / T;
E = h_bar ^ 2 * k .^ 2 / (2 * m0) + Vnumber(n_size + 1);


% 利用微扰论进行修正，此时为扩展布里渊区图景
for k = -n_size : acc : n_size
    N = round(k);
    n = abs(N);
    if(n ~= 0 && abs(k - N) < p && Vnumber(N + n_size + 1) ~= 0 && k <= 0)
        %对布里渊区边界附近做简并修正
        Ek = E(round((k + n_size) / acc) + 1);
        Ekk = E(round((k + 2 * n + n_size) / acc) + 1);
        E1 = 0.5 * (Ek + Ekk + ((Ek - Ekk) ^ 2 + 4 *...
            Vnumber(N + n_size + 1) ^ 2) ^ 0.5);
        E2 = 0.5 * (Ek + Ekk - ((Ek - Ekk) ^ 2 + 4 *...
            Vnumber(N + n_size + 1) ^ 2) ^ 0.5);
        if(Ek > Ekk)    %能级劈裂
            E(round((k + n_size) / acc) + 1) = E1;
            E(round((k + 2 * n + n_size) / acc) + 1) = E2;
        else
            E(round((k + n_size) / acc) + 1)= E2;
            E(round((k + 2 * n + n_size) / acc) + 1) = E1;
        end       
    else
        if(n == 0 || (abs(k - N) >= p) || Vnumber(N + n_size + 1) == 0)
            %需要做非简并修正，只采用二阶修正
            Ek = E(round((k + n_size) / acc) + 1);
            for i = -n_size : n_size
                if(i ~= 0)
                    Ek = Ek + Vnumber(i + n_size + 1) ^ 2 / (h_bar ^ 2 /...
                    (2 * m0) * (pi / T) ^ 2 * (k ^ 2 - (k + 2 * i) ^ 2));
                end
            end
            E(round((k + n_size) / acc) + 1) = Ek;
        end
    end
end

% 通过平移化为简约布里渊区图景
l = length(-1 : acc : 1);
e = zeros(n_size - 1, l);
for k = 1 : n_size - 1
    if mod(k, 2) == 1
        e(k, :) = [E(round((-k + n_size) / acc) + 1 :...
            round((-k + n_size + 1) / acc)), ...
            E(round((k + n_size - 1) / acc) + 1 :...
            round((k + n_size) / acc)),...
            E(round((-k + n_size) / acc) + 1)];
    end
    if mod(k, 2) == 0
        e(k, :) = [E(round((k + n_size - 1) / acc) + 1 :...
            round((k + n_size) / acc)), ...
            E(round((-k + n_size) / acc) + 1 :...
            round((-k + n_size + 1) / acc)),...
            E(round((k + n_size - 1) / acc) + 1)];
    end
end
