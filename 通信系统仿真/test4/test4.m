%% 4-1
% 定义信号参数
A = 1; % 振幅
a = 0.5; % 指数因子
omega_0 = 2 * pi; % 基频
phi = pi / 4; % 相位

% 定义时间范围
t = linspace(0, 10, 1000); % 从0到10，生成1000个点

% 定义信号
y = A * exp(a * t) .* cos(omega_0 * t + phi);

% 绘制波形
figure;
plot(t, y);
title('信号 Ae^{at} \\cos(\\omega_0 t + \\phi) 的波形'); % 修复了标题的特殊字符
xlabel('时间 (s)');
ylabel('幅度');
grid on;


%% 4-2
s = tf('s'); % 定义拉普拉斯变量

% 定义传递函数
H = 10 / (s^2 + 2 * s + 100); % 系统的传递函数
% 输入信号的拉普拉斯变换
F_s = 2 * pi / (s^2 + (2 * pi)^2); % sin(2 * pi * t) 的拉普拉斯变换
% 输出信号的拉普拉斯变换
Y_s = H * F_s; % 传递函数与输入信号拉普拉斯变换的乘积

% 计算零状态响应的时间域
t = 0:0.01:10; % 定义时间范围
y_t = impulse(Y_s, t); % 使用冲激响应函数求逆拉普拉斯变换

% 绘制零状态响应
figure;
plot(t, y_t);
title('系统零状态响应');
xlabel('时间 (s)');
ylabel('响应');
grid on;

%% 4-3
s = tf('s'); % 定义拉普拉斯变量

% 定义系统的传递函数
H = 10 / (s^2 + 2 * s + 100); % 系统的传递函数
% 计算冲激响应
t = 0:0.01:10; % 时间范围
y_t = impulse(H, t); % 计算冲激响应

% 绘制冲激响应波形
figure;
plot(t, y_t);
title('系统冲激响应');
xlabel('时间 (s)');
ylabel('响应幅度');
grid on; % 显示网格线

%% 4-4
% 定义 Z 变量
z = tf('z', 1); % 定义 Z 变换变量，采样时间为1

% 定义系统传递函数
H = 10 / (1 + 3/z + 2/(z^2)); % 系统的传递函数

% 计算单位脉冲响应
N = 10; % 定义计算的长度
y_impulse = impulse(H, N); % 计算单位脉冲响应

% 重新定义 X 轴
x_values = 0:(N-1); % X 轴数据


% 确保两者的长度相等
max_length = max(length(x_values), length(y_impulse)); % 找到最大长度
x_values = [x_values, (length(x_values):(max_length-1))]; % 填充 x_values
y_impulse = [y_impulse; zeros(max_length - length(y_impulse), 1)]; % 填充 y_impulse

% 绘制单位脉冲响应
figure;
stem(x_values, y_impulse); % 确保 X 和 Y 长度一致
title('系统单位脉冲响应');
xlabel('k');
ylabel('响应幅度');
grid on;



%% 4-5
% 定义时间范围
t=-6*pi:0.0001:6*pi;

% 定义周期为 2 的三角波，幅度为 -1 到 1
tri_wave=-sawtooth(pi*(t-0.5),0.5);

% 绘制三角波
figure;
plot(t, tri_wave);
title('周期为 2 的三角波信号');
xlabel('时间 (s)');
ylabel('幅度');
grid on;

% 计算 FFT
N = length(t); % 信号的长度
Fs = N / (t(end) - t(1)); % 采样频率
tri_wave_fft = fft(tri_wave); % 计算 FFT

% 计算频率
f = (0:(N/2 - 1)) * (Fs / N); % 频率范围，只看正频率部分

% 计算 FFT 的幅度
fft_magnitude = abs(tri_wave_fft / N); % 归一化 FFT 幅度
fft_magnitude = fft_magnitude(1:(N/2)); % 只看正频率部分

% 绘制频谱
figure;
plot(f, fft_magnitude);
title('三角波信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅度');
grid on; % 显示网格

%% 4-6
% 定义时间范围和采样点数
T = 2; % 信号周期
Fs = 100; % 采样频率
t = 0:1/Fs:4*T; % 时间范围，生成足够多的点

% 生成三角波信号，周期为 2，幅度 ±1
tri_wave = 2 * sawtooth(2 * pi / T * (t-1), 0.5) - 1; % 幅度调整到 -1 到 1

% 绘制三角波信号
figure;
plot(t, tri_wave);
title('三角波信号');
xlabel('时间 (s)');
ylabel('幅度');
grid on;

% 计算 FFT
N = length(tri_wave); % 信号的长度
fft_tri_wave = fft(tri_wave); % 应用 FFT

% 计算频率轴
f = (0:N-1) * (Fs / N); % 频率范围

% 计算幅度谱
fft_magnitude = abs(fft_tri_wave / N); % 归一化 FFT 幅度

% 绘制频谱图，只显示正频率部分
figure;
plot(f(1:floor(N/2)), fft_magnitude(1:floor(N/2)));
title('三角波信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅度');
grid on;

%% 4-7
% 定义时间范围
t = linspace(0, 10, 1000); % 时间范围，生成 1000 个点

% 定义周期方波信号
Vin = square(2*pi*t); % 周期方波信号

% 定义RC电路参数
R = 1; % 电阻（Ω）
C = 1; % 电容（F）
RC = R * C; % RC 时间常数

% 计算系统的频率响应
s = tf('s'); % 创建传递函数变量
H = 1 / (RC*s + 1); % RC系统的传递函数

% 计算周期方波信号通过RC系统的响应
Vout = lsim(H, Vin, t); % 使用 lsim 函数计算响应

% 绘制输入输出信号
figure;
plot(t, Vin, 'b', t, Vout, 'r');
title('周期方波通过RC系统的响应');
xlabel('时间 (s)');
ylabel('电压 (V)');
legend('输入信号', '输出信号');
grid on;

%% 4-8
% 定义参数
alpha = 0.5; % 可调整 alpha 值

% 定义传递函数的分子和分母
b = 1; % 分子系数
a = [1, -alpha]; % 分母系数

% 使用 freqz 计算频率响应
[H_freq, w] = freqz(b, a, 512); % 512 是计算点数

% 绘制幅度响应曲线
figure;
plot(w, abs(H_freq)); % 绘制幅度
title('H(e^{j\Omega}) 的幅度响应曲线');
xlabel('频率 (rad/s)');
ylabel('幅度');
grid on;

