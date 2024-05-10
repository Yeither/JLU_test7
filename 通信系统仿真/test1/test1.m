%% 1-2
A=[7  1  5;2  5  6;3  1  5],B=[1  1  1; 2  2  2; 3  3  3]
% 获取第2行第3列的元素
result1 = A(2,3) % 输出：6

% 获取矩阵A的第2列
result2 = A(:,2) % 输出：[1; 5; 1]

% 获取矩阵A的第3行
result3 = A(3,:) % 输出：[3 1 5]

% 获取矩阵A的第1和第3列
result4 = A(:,1:2:3) % 输出：[7 5; 2 6; 3 5]

% 点乘操作，A的第3列和B的第2列
result5 = A(:,3) .* B(:,2) % 输出：[5; 12; 15]

% 矩阵乘法，A的第3列乘以B的第2行
result6 = A(:,3) * B(2,:) % 输出：[10 10 10; 12 12 12; 15 15 15]

% 矩阵乘法，A与B
result7 = A * B % 输出：[22 22 22; 35 35 35; 23 23 23]

% 点乘操作，A与B
result8 = A .* B % 输出：[7 1 5; 4 10 12; 9 3 15]

% 矩阵幂运算，A与自身相乘
result9 = A^2 % 输出：[58 38 70; 38 43 67; 29 23 48]

% 元素幂运算，A的每个元素平方
result10 = A.^2 % 输出：[49 1 25; 4 25 36; 9 1 25]

% 矩阵除法，B除以A
result11 = B / A % 输出：[0.2361 0.2778 -0.1667; 0.4444 0.3889 0.1667; 0.3194 0.3333 -0.0556]

% 元素除法，B的每个元素除以A的对应元素
result12 = B ./ A % 输出：[0.1429 1 0.2; 1 0.4 0.3333; 1 3 0.2]

%% 1-3
C = 1:2:20; % 生成 [1 3 5 7 9 11 13 15 17 19]

%% 1-4
% 创建一些变量
A = [1, 2, 3];
B = magic(3);
C = 'Hello, World!';

% 查看所有变量信息
whos % 显示当前工作区中所有变量的详细信息

% 假设我们不需要变量 C
clear C % 删除变量 C

% 查看剩余变量信息
whos % 显示剩余的变量信息

%% 1-5
% 定义变量 x 和 y
x = 2;
y = 4;

% 计算表达式
z = x^2 + exp(x + y) - y * log(x) - 3;

% 显示结果
disp("z 的值是: " + z);

%% 1-6-1
x = -2:0.01:2; % 定义x的范围
y = x.^2 .* sin(x.^2 - x - 2); % 计算y(x)

% 使用plot绘制图像
figure;
plot(x, y);
title('y(x) = x^2 * sin(x^2 - x - 2) (Using plot)');
xlabel('x');
ylabel('y');

% 使用fplot绘制图像
figure;
fplot(@(x) x.^2 .* sin(x.^2 - x - 2), [-2, 2]);
title('y(x) = x^2 * sin(x^2 - x - 2) (Using fplot)');
xlabel('x');
ylabel('y');

%% 1-6-2
% 定义参数
t = linspace(0, 2*pi, 100);

% 椭圆参数方程
x = 2 * cos(t);
y = 4 * sin(t);

% 绘制椭圆
figure;
plot(x, y);
title('x^2/4 + y^2/16 = 1');
xlabel('x');
ylabel('y');
axis equal; % 确保轴的比例相同

%% 1-6-3
% 定义x和y的范围
x = -3:0.1:3;
y = -2:0.1:2;

% 创建网格
[X, Y] = meshgrid(x, y);

% 定义函数z
Z = (X.^2 - 2 * X) .* exp(-(X.^2 + Y.^2 + X .* Y));

% 绘制3D网格图
figure;
mesh(X, Y, Z);
title('3D Mesh Grid of z=(x^2-2x) * exp(-(x^2 + y^2 + xy))');
xlabel('x');
ylabel('y');
zlabel('z');

% 绘制二维等高线图
figure;
contour(X, Y, Z, 20); % 20个等高线
title('2D Contour of z=(x^2-2x) * exp(-(x^2 + y^2 + xy))');
xlabel('x');
ylabel('y');

% 绘制三维等高线图
figure;
contour3(X, Y, Z, 20); % 20个等高线
title('3D Contour of z=(x^2-2x) * exp(-(x^2 + y^2 + xy))');
xlabel('x');
ylabel('y');
zlabel('z');

%% 1-6-4
x = linspace(0, 2*pi, 100);

% 定义不同的函数
y1 = cos(x);
y2 = sin(x - pi/2);
y3 = x.^2 .* cos(x - pi);
y4 = exp(sin(x));

% 创建子图
figure;

% 子图1
subplot(2, 2, 1); % 2x2 网格，位置 1
plot(x, y1);
title('y1 = cos(x)');

% 子图2
subplot(2, 2, 2); % 2x2 网格，位置 2
plot(x, y2);
title('y2 = sin(x - pi/2)');

% 子图3
subplot(2, 2, 3); % 2x2 网格，位置 3
plot(x, y3);
title('y3 = x^2 * cos(x - pi)');

% 子图4
subplot(2, 2, 4); % 2x2 网格，位置 4
plot(x, y4);
title('y4 = exp(sin(x))');

%% 1-7
% 定义x的范围
x = linspace(0, 2, 100); % 从0到2，100个点

% 计算y的不同函数
y1 = sqrt(x); % y = √x
y2 = x.^2;   % y = x^2
y3 = nthroot(x, 3); % y = ∛x
y4 = x.^3;   % y = x^3
y5 = x;      % y = x

% 绘制图像
figure;
hold on; % 确保在同一图中叠加绘制

plot(x, y1, 'b-', 'LineWidth', 1.5); % 蓝色实线
plot(x, y2, 'r--', 'LineWidth', 1.5); % 红色虚线
plot(x, y3, 'g-.', 'LineWidth', 1.5); % 绿色点划线
plot(x, y4, 'm:', 'LineWidth', 1.5); % 紫色点点线
plot(x, y5, 'k-', 'LineWidth', 1.5); % 黑色实线

% 设置图表标题和轴标签
title('不同函数的图像');
xlabel('x');
ylabel('y');

% 添加图例
legend('y = √x', 'y = x^2', 'y = ∛x', 'y = x^3', 'y = x');

% 显示网格
grid on;

hold off; % 结束绘制
