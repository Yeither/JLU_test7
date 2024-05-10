%% 7-1-1
% 定义传递函数
% 简化分子多项式
num = [4, 8, 4, 96, 144, 144];  % 降低分子阶数
den = [1, 3, 5, 6, 2, 5];

G = tf(num, den);

% 零极点图
figure;
pzmap(G);
title('零极点图');

% 频率响应 - Bode 图
figure;
bode(G);
title('Bode 图');

% 时间响应 - 阶跃响应
figure;
step(G);
title('阶跃响应');

%% 7-1-2
% 定义传递函数
num = [4, -2];
den = [1, 0, 2, 5];
G = tf(num, den);

% 零极点图
figure;
pzmap(G);
title('零极点图');

% 频率响应 - Bode 图
figure;
bode(G);
title('Bode 图');

% 时间响应 - 阶跃响应
figure;
step(G);
title('阶跃响应');

%% 7-1-3
% 定义传递函数
num = [1];
den = [2, 5, 2];
G = tf(num, den);

% 零极点图
figure;
pzmap(G);
title('零极点图');

% 频率响应 - Bode 图
figure;
bode(G);
title('Bode 图');

% 时间响应 - 阶跃响应
figure;
step(G);
title('阶跃响应');

%% 7-1-4
% 定义状态空间模型
A = [0, 1; -1, -2];
B = [0; 1];
C = [0, 1];
D = [0];

% 转换为传递函数
[num, den] = ss2tf(A, B, C, D);
G = tf(num, den);

% 零极点图
figure;
pzmap(G);
title('零极点图');

% 频率响应 - Bode 图
figure;
bode(G);
title('Bode 图');

% 时间响应 - 阶跃响应
figure;
step(G);
title('阶跃响应');

%% 7-2
% 定义开环传递函数
num = [1];  % k 在根轨迹图上是增益，通常不指定
den = [0.5, 1, 0, 4, 1];
G = tf(num, den);

% 绘制根轨迹图
figure;
rlocus(G);
title('根轨迹图');

%% 7-3
% 定义前向通路传递函数
G = tf(10, [1, 0, 1]);

% 定义反馈通路传递函数
H_f = tf([0.3, 1], 1);

% 计算闭环传递函数
T = feedback(G, H_f);  % 负反馈系统

% 绘制单位阶跃响应曲线
figure;
[response, time] = step(T);  % 获取阶跃响应
plot(time, response);
title('单位阶跃响应');
xlabel('时间 (秒)');
ylabel('响应');
grid on;

% 计算上升时间、峰值时间、超调量
S = stepinfo(T);

% 上升时间、峰值时间、超调量
disp('上升时间: ');
disp(S.RiseTime);

disp('峰值时间: ');
disp(S.PeakTime);

disp('超调量: ');
disp(S.Overshoot);

% 计算延迟时间
threshold = 0.1 * response(end);  % 阈值设为稳态的10%
delay_time = time(find(response >= threshold, 1));  % 首次超过阈值的时间

disp('延迟时间: ');
disp(delay_time);
