%% 5-1
% Define parameters
n = 0:15;  % Time indices

% Parameters for Gaussian sequence
p_default = 8;
q_values = [2, 4, 8];
p_values = [8, 13, 14];
a = 0.1;  % For decaying sine
f = 0.1;  % Frequency

% Generate Gaussian sequences with varying q
for q = q_values
    x_a = exp(-((n - p_default).^2) / q);
    plot_time_freq(x_a, ['Gaussian sequence with q=', num2str(q)]);
end

% Generate Gaussian sequences with varying p
for p = p_values
    x_a = exp(-((n - p).^2) / q_values(end)); % Use the last q as default
    plot_time_freq(x_a, ['Gaussian sequence with p=', num2str(p)]);
end

%% 5-2
% 参数定义
n = 0:15; % 定义时间范围
a = 0.2; % 衰减因子
f1 = 0.1; % 正弦频率 1
f2 = 0.5; % 正弦频率 2

% 衰减正弦序列
x_b1 = exp(-a * n) .* sin(2 * pi * f1 * n); % 频率为f1
x_b2 = exp(-a * n) .* sin(2 * pi * f2 * n); % 频率为f2

% 绘制时域信号
figure;
subplot(2, 1, 1);
stem(n, x_b1);
title('时域信号 - 衰减正弦序列 (f = 0.1)');
xlabel('n');
ylabel('x_b(n)');
grid on;

subplot(2, 1, 2);
stem(n, x_b2);
title('时域信号 - 衰减正弦序列 (f = 0.5)');
xlabel('n');
ylabel('x_b(n)');
grid on;

% 计算并绘制频谱
X_b1 = fft(x_b1, 256); % 零填充至256点FFT
X_b2 = fft(x_b2, 256); % 零填充至256点FFT

f = (0:127) / 256; % 频率范围

figure;
subplot(2, 1, 1);
plot(f, abs(X_b1(1:128))); % 绘制正频率范围
title('频域信号 - 衰减正弦序列 (f = 0.1)');
xlabel('频率');
ylabel('|X_b(f)|');
grid on;

subplot(2, 1, 2);
plot(f, abs(X_b2(1:128))); % 绘制正频率范围
title('频域信号 - 衰减正弦序列 (f = 0.5)');
xlabel('频率');
ylabel('|X_b(f)|');
grid on;

%% 5-3
% 定义序列范围
n = 0:15; % 序列长度

% 三角波序列
x_c = [n(1:4) + 1, 8 - n(5:8)]; % 三角波序列
x_c = [x_c, zeros(1, 8)]; % 在末尾补零

% 反三角波序列
x_d = [4 - n(1:4), n(5:8) - 3]; % 反三角波序列
x_d = [x_d, zeros(1, 8)]; % 在末尾补零

% 绘制时域信号
figure;
subplot(2, 1, 1);
stem(0:15, x_c(1:16));
title('三角波序列 - 时域');
xlabel('n');
ylabel('x_c(n)');
grid on;

subplot(2, 1, 2);
stem(0:15, x_d(1:16));
title('反三角波序列 - 时域');
xlabel('n');
ylabel('x_d(n)');
grid on;

% 使用256点FFT进行频域分析
N = 256; % FFT的点数
X_c = fft(x_c, N); % 计算三角波序列的FFT
X_d = fft(x_d, N); % 计算反三角波序列的FFT

% 绘制频域信号
f = (0:(N/2 - 1)) * (1 / N); % 正频率范围

figure;
subplot(2, 1, 1);
plot(f, abs(X_c(1:N/2)));
title('三角波序列 - 幅频特性');
xlabel('频率 (单位: 周期/采样点数)');
ylabel('|X_c(f)|');
grid on;

subplot(2, 1, 2);
plot(f, abs(X_d(1:N/2)));
title('反三角波序列 - 幅频特性');
xlabel('频率 (单位: 周期/采样点数)');
ylabel('|X_d(f)|');
grid on;


%% 5-4
% 参数定义
N1 = 16; % 第一个样本数
N2 = 128; % 第二个样本数
f_base = 0.125; % 基础频率

% Δf 分别为 1/16 和 1/64
delta_f1 = 1/16; % 第一个 Δf
delta_f2 = 1/64; % 第二个 Δf

% 生成信号
n1 = 0:(N1 - 1); % 时间索引
x1 = sin(2 * pi * f_base * n1) + cos(2 * pi * (f_base + delta_f1) * n1);

n2 = 0:(N2 - 1); % 时间索引
x2 = sin(2 * pi * f_base * n2) + cos(2 * pi * (f_base + delta_f2) * n2);

% 计算频谱
X1 = fft(x1); % 对第一个信号进行FFT
X2 = fft(x2); % 对第二个信号进行FFT

% 频率轴
f1 = (0:(N1 - 1)) / N1; % 对应N1的频率轴
f2 = (0:(N2 - 1)) / N2; % 对应N2的频率轴

% 绘制频谱
figure;
subplot(2, 1, 1);
plot(f1, abs(X1)); % 绘制第一个信号的频谱
title('频谱 - N = 16, Δf = 1/16');
xlabel('归一化频率');
ylabel('|X(f)|');
grid on;

subplot(2, 1, 2);
plot(f2, abs(X2)); % 绘制第二个信号的频谱
title('频谱 - N = 128, Δf = 1/64');
xlabel('归一化频率');
ylabel('|X(f)|');
grid on;

%% 5-5
% 定义参数
p = 8; % Gaussian序列的均值
q = 2; % Gaussian序列的方差
a = 0.1; % 衰减正弦序列的衰减因子
f = 0.0625; % 衰减正弦序列的频率

% 定义序列长度
n = 0:15; % 时间范围

% 定义Gaussian序列
x_a = exp(-((n - p).^2) / q); % Gaussian序列

% 定义衰减正弦序列
x_b = exp(-a * n) .* sin(2 * pi * f * n); % 衰减正弦序列

% 圆周卷积（16点FFT）
X_a = fft(x_a, 16); % 对x_a进行16点FFT
X_b = fft(x_b, 16); % 对x_b进行16点FFT

% FFT的乘积实现圆周卷积
circ_conv_result = ifft(X_a .* X_b); % 圆周卷积

% 线性卷积（需要零填充以避免循环效应）
X_a_padded = fft(x_a, 31); % 零填充至31点，保证线性卷积长度
X_b_padded = fft(x_b, 31); % 零填充至31点

% FFT的乘积实现线性卷积
linear_conv_result = ifft(X_a_padded .* X_b_padded); % 线性卷积

% 绘制结果
figure;
subplot(2, 1, 1);
stem(0:15, circ_conv_result);
title('圆周卷积 - x_a(n) 和 x_b(n)');
xlabel('n');
ylabel('结果');
grid on;

subplot(2, 1, 2);
stem(0:30, linear_conv_result);
title('线性卷积 - x_a(n) 和 x_b(n)');
xlabel('n');
ylabel('结果');
grid on;


%% 5-6

% 生成高斯序列 x_a
p = 7.5; % 中心位置
q = 2; % 标准差
n = 0:15; % 序列索引
x_a = exp(-((n - p).^2) / q);

% 生成衰减正弦序列 x_b
a_n = 0.3; % 衰减因子
f_n = 0.2; % 正弦频率
x_b = exp(-a_n * n) .* sin(2 * pi * f_n * n);

% 生成三角波序列 x_c
x_c = [n(1:4) + 1, 8 - n(5:8)];

% 生成反三角波序列 x_d
x_d = [4 - n(1:4), n(5:8) - 3];

% 产生 512 点随机序列 x_rand
x_rand = randn(1, 512);

% 对 x_rand 与 x_c 进行卷积
y_conv_c = conv(x_rand, x_c);

% 对 x_rand 与 x_d 进行卷积
y_conv_d = conv(x_rand, x_d);

% 重叠相加法
segment_size = 64; % 512/8 = 64
y_overlap_add = zeros(1, 512 + length(x_c) - 1); % 创建足够大的结果数组

% 遍历 8 段
for i = 0:7
    start_index = i * (segment_size - (length(x_c) - 1)); % 注意重叠部分
    segment = x_rand(start_index + 1 : start_index + segment_size); % 获取分段
    conv_segment = conv(segment, x_c); % 进行卷积
    % 确保索引范围在允许范围内
    end_index = start_index + length(conv_segment);
    if end_index <= length(y_overlap_add)
        % 重叠相加
        y_overlap_add(start_index + 1 : end_index) = ...
            y_overlap_add(start_index + 1 : end_index) + conv_segment;
    end
end

% 绘制重叠相加法的结果
figure;
plot(y_overlap_add);
title('重叠相加法的结果');
xlabel('样本');
ylabel('幅度');



%% 5-7
% 定义滤波器规格
f_p = 0.2e3; % 通带频率 (0.2 kHz)
f_r = 0.3e3; % 阻带频率 (0.3 kHz)
A_p = 1; % 通带衰减 (1 dB)
A_r = 25; % 阻带衰减 (25 dB)
T = 1e-3; % 采样时间 (1 ms)

% 转换频率至归一化角频率 (rad/s)
omega_p = 2 * pi * f_p; % 通带角频率
omega_r = 2 * pi * f_r; % 阻带角频率

% 计算Butterworth滤波器的阶数和截止频率
[N, Wn] = buttord(omega_p, omega_r, A_p, A_r, 's'); % 计算模拟滤波器的阶数和截止频率
[b, a] = butter(N, Wn, 's'); % 模拟低通Butterworth滤波器的传递函数

% 脉冲响应不变法
[b_ii, a_ii] = impinvar(b, a, 1 / T); % 通过脉冲响应不变法转换至数字滤波器

% 双线性变换法
[b_bl, a_bl] = bilinear(b, a, 1 / T); % 通过双线性变换转换至数字滤波器

% 绘制幅频特性
fvtool(b_ii, a_ii, b_bl, a_bl); % 使用滤波器可视化工具比较两个滤波器
legend('脉冲响应不变法', '双线性变换法');

%% 5-1
% Helper function to plot time domain and frequency domain
function plot_time_freq(seq, title_str)
    % Plot the time-domain sequence
    figure;
    subplot(2, 1, 1);
    stem(seq, 'filled');
    title(['Time-domain: ', title_str]);
    xlabel('n'); ylabel('Amplitude');
    
    % Plot the frequency-domain sequence (magnitude spectrum)
    subplot(2, 1, 2);
    mag_spectrum = abs(fft(seq, 256)); % FFT with zero-padding
    plot(linspace(0, 1, 256), mag_spectrum);
    title(['Frequency-domain: ', title_str]);
    xlabel('Normalized Frequency'); ylabel('Magnitude');
end
