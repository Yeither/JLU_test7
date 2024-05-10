%% 2-1
% 输入 n 的值
n = input('请输入 n 的值: ');

% 计算等差数列的和
% 首项 a = 1，公差 d = 2，项数 m = n + 1
% 求和公式: S = m/2 * (2 * a + (m - 1) * d)
m = n + 1; % 项数
a = 1; % 首项
d = 2; % 公差

S = m / 2 * (2 * a + (m - 1) * d); % 求和

% 输出结果
disp(['1 + 3 + 5 + ... + (2 * ', num2str(n), ' + 1) = ', num2str(S)]);

%% 2-2
% 计算 f(-3)
result1 = ff(-3);
disp(['f(-3) = ', num2str(result1)]);

% 计算 f(√2)
result2 = ff(sqrt(2));
disp(['f(√2) = ', num2str(result2)]);

% 计算 f(∞)
result3 = ff(Inf);
disp(['f(∞) = ', num2str(result3)]);


%% 2-3
% 定义矩阵 A
A = [1 -6 3 2; 3 -5 4 0; -1 -11 2 4];

% 求矩阵 A 的秩
rank_A = rank(A); % 使用 MATLAB 的 rank 函数

% 输出结果
disp(['Rank(A) = ', num2str(rank_A)]);

%% 2-4-1
% 输入旋转角度
t = input('请输入旋转角度 t（弧度）: ');

% 构建旋转矩阵
A = [cos(t), -sin(t); sin(t), cos(t)];

% 计算五次幂 A^5
theta = 5 * t; % 旋转角度乘以 5
A_5 = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% 输出结果
disp('A^5 =');
disp(A_5);

%% 2-4-2
% 定义矩阵 B
B = [
    1, 2, 1, 0;
    6, 2, 4, 1;
    0, 2, 1, 0;
    3, 1, 4, 1
];

% 求逆矩阵 B^(-1)
B_inv = inv(B); % 使用 inv 函数

% 输出结果
disp('B^(-1) =');
disp(B_inv);

%% 2-5
% 输入变量
x = 1;
n = 100;

% 计算 e^x 的近似值
e_approx = approx_exp(x, n);

% 输出结果
disp(['e^', num2str(x), ' 的近似值是：', num2str(e_approx)]);

%% 2-6-1
% 多项式 f(x) 的系数
f_coeff = [3, -1, 2, 1, 1, 3];

% 求 f(x) 的根
f_roots = roots(f_coeff);

% 输出根
disp('f(x) 的根是:');
disp(f_roots);

%% 2-6-2
% 多项式 g(x) 的系数
g_coeff = [1/3, 1, -3, -1];

% 定义多项式 g(x)
g = @(x) polyval(g_coeff, x);

% 求 g(x) 在区间 [-1, 2] 上的最小值
[x_min, g_min] = fminbnd(g, -1, 2);

% 输出最小值
disp(['g(x) 在 [-1, 2] 上的最小值是: ', num2str(g_min), ', 在 x = ', num2str(x_min), ' 处。']);

%% 2-6-3
% 定义多项式 f(x) 和 g(x) 的系数
f_coeff = [3, -1, 2, 1, 1, 3];
g_coeff = [1/3, 1, -3, -1];

% 确保系数长度一致，通过填充 0
len_f = length(f_coeff);
len_g = length(g_coeff);
if len_f > len_g
    g_coeff = [zeros(1, len_f - len_g), g_coeff];
elseif len_g > len_f
    f_coeff = [zeros(1, len_g - len_f), f_coeff];
end

% 求 f(x) + g(x)
fg_sum = f_coeff + g_coeff; % 多项式加法

% 输出结果
disp('f(x) + g(x) 的系数是:');
disp(fg_sum);
% 计算 f(x) ⋅ g(x)
fg_prod = conv(f_coeff, g_coeff);

% 输出结果
disp('f(x) ⋅ g(x) 的系数是:');
disp(fg_prod);
% 计算 f(x) / g(x)
% 确保 g_coeff 的首项非零
if g_coeff(1) == 0
    % 找到第一个非零系数
    idx = find(g_coeff ~= 0, 1);
    % 去掉前面的零
    g_coeff = g_coeff(idx:end);
end

% 使用 deconv 进行多项式除法
[fg_quotient, fg_remainder] = deconv(f_coeff, g_coeff);

% 输出结果
disp('f(x) / g(x) 的商是:');
disp(fg_quotient);

disp('f(x) / g(x) 的余数是:');
disp(fg_remainder);

% 输出结果
disp('f(x) / g(x) 的商是:');
disp(fg_quotient);

disp('f(x) / g(x) 的余数是:');
disp(fg_remainder);


%% 2-6-4
% 求 f(x) 的导数
f_derivative = polyder(f_coeff);

% 输出 f(x) 的导数
disp('f(x) 的导数是:');
disp(f_derivative);

%% 2-7
% 计算64个格子的麦子总数
total_grains = 2^64 - 1;

% 每袋小麦约含 1.4 × 10^8 粒
grains_per_bag = 1.4e8;

% 计算需要多少袋小麦
bags_needed = total_grains / grains_per_bag;

% 输出结果
disp(['总共需要小麦数量: ', num2str(total_grains)]);
disp(['需要多少袋小麦: ', num2str(bags_needed)]);

%% 2-8
% 程序计算满足条件的鸡翁、鸡母、鸡雏的数量
% 鸡翁5钱一只，鸡母3钱一只，鸡雏三只1钱

% 为了满足"百钱买百鸡"，我们可以分别枚举鸡翁和鸡母的数量
% 鸡翁数量最大不能超过20，因为20只公鸡要花100钱
for rooster = 0:20
    % 鸡母数量最大不能超过33，因为33只母鸡刚好100钱
    for hen = 0:33
        % 根据鸡翁和鸡母数量计算鸡雏数量
        chick = 100 - rooster - hen; % 剩余的鸡数
        % 根据鸡雏数量计算所需的费用
        total_cost = rooster * 5 + hen * 3 + (chick / 3);
        
        % 检查费用是否正好100钱，并且鸡雏数量是3的倍数
        if total_cost == 100 && mod(chick, 3) == 0
            fprintf('鸡翁: %d, 鸡母: %d, 鸡雏: %d\n', rooster, hen, chick);
        end
    end
end


%% 2-2
function y = ff(x)
% 分段函数 f(x)={■(x&0≤x<1@2-x&1≤x≤2@0&其它)┤

if x >= 0 && x < 1
    y = x;
elseif x >= 1 && x <= 2
    y = 2 - x;
else
    y = 0; % 对于其它情况，返回 0
end
end

%% 2-5
function result = approx_exp(x, n)
% 计算 e^x 的近似值，使用泰勒级数展开
% 输入 x 为自变量，n 为展开的最大阶数

% 初始化结果
result = 1; % 第一项是 1
term = 1;   % 用于计算每一项的中间变量

% 循环计算每一项，并累加
for k = 1:n
    term = term * x / k; % 计算第 k 项
    result = result + term; % 累加到结果
end
end
