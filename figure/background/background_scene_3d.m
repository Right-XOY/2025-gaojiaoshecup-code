%% background_scene_3d.m
% 2025 高教社杯 A 题「烟幕干扰弹的投放策略」背景介绍三维示意图
% 目的：直观展示 来袭导弹(M1~M3)、无人机(FY1~FY5)、假目标(原点)、
%       真目标(圆柱) 的空间关系，用于论文"问题背景/场景描述"部分。
%
% 布局：分两张图输出
%       (a) 全局态势 —— 三枚导弹与五架无人机的初始位置及来袭方向；
%       (b) 目标区放大 —— 假目标原点 + 真目标圆柱(半径7m、高10m)，
%           体现真目标相对假目标在 +Y 方向 200 m 的偏移。
%
% 输出：background_scene_overview.png（全局态势）、
%       background_scene_target.png（目标区放大），均 300 dpi

clear; clc; close all;

%% ================== 1. 场景常量与坐标 ==================
% 来袭导弹初始位置（3x3，每列一枚：M1/M2/M3）
M0 = [20000 19000 18000;
         0   600  -600;
      2000  2100  1900];
m_name = {'M1','M2','M3'};

% 无人机初始位置（3x5，每列一架：FY1~FY5）
U0 = [17800 12000  6000 11000 13000;
         0  1400 -3000  2000 -2000;
      1800  1400   700  1800  1300];
u_name = {'FY1','FY2','FY3','FY4','FY5'};

% 假目标（导弹瞄准点）与真目标（圆柱形固定目标）
fake   = [0; 0; 0];        % 假目标 = 原点，导弹飞行方向直指此处
true_c = [0; 200; 0];      % 真目标下底面圆心
Rc = 7;  Hc = 10;          % 真目标圆柱半径、高

% 配色
c_m   = [0.85 0.20 0.20];  % 导弹：红
c_u   = [0.15 0.40 0.75];  % 无人机：蓝
c_fake= [0.15 0.40 0.75];  % 假目标：蓝（同无人机）
c_true= [0.15 0.40 0.75];  % 真目标：蓝（同无人机）
fs    = 10;                % 字号

%% ================== 2. 图 (a)：全局态势 ==================
figure('Color','w','Position',[60 60 720 500]);
hold on; grid on; view([-38, 26]);
% 不采用 axis equal：X 轴跨度 20 km、Y 5 km、Z 2.1 km，等比会使图被压成极扁的条带。
% 用 daspect 压缩 X、适度拉伸 Z，使整体展宽/高度更均衡（数值可按需微调）。
daspect([6 2 1]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

% 导弹：红色三角 + 指向假目标的来袭方向（虚线）
% 各导弹编号在 Y 方向的微调（M3 往 -Y 移开，避免与 M1/M2 或方向线重叠）
m_dy = [0 0 -250];
for i = 1:3
    plot3([M0(1,i) fake(1)], [M0(2,i) fake(2)], [M0(3,i) fake(3)], ...
          '--', 'Color', c_m, 'LineWidth', 1.1);
    scatter3(M0(1,i), M0(2,i), M0(3,i), 95, c_m, '^', 'filled');
    text(M0(1,i), M0(2,i)+m_dy(i), M0(3,i)+200, m_name{i}, ...
         'Color', c_m, 'FontSize', fs, 'FontWeight', 'bold');
end

% 无人机：蓝色圆点
% 各无人机编号在 X/Y 方向的微调（FY1 往 +Y、-X 移开）
u_dx = [-300 0 0 0 0];
u_dy = [ 300 0 0 0 0];
for i = 1:5
    scatter3(U0(1,i), U0(2,i), U0(3,i), 70, c_u, 'o', 'filled');
    text(U0(1,i)+u_dx(i), U0(2,i)+u_dy(i), U0(3,i)+200, u_name{i}, ...
         'Color', c_u, 'FontSize', fs);
end

% 假目标（原点）：蓝色圆点
scatter3(0, 0, 0, 120, c_fake, 'o', 'filled');
text(3000, 0, -220, '假目标', 'Color', c_fake, 'FontSize', fs, ...
     'HorizontalAlignment', 'center');

set(gca, 'FontName', 'Times New Roman', 'FontSize', fs, 'LineWidth', 0.8, ...
    'XTick', 0:5000:20000, 'XTickLabel', {'0','5000','10000','15000','20000'}, ...
    'YTick', -3000:1000:2000, 'YTickLabel', {'-3000','-2000','-1000','0','1000','2000'}, ...
    'ZTick', 0:500:2000, 'ZTickLabel', {'0','500','1000','1500','2000'});
box on;

print(gcf, 'background_scene_overview', '-dpng', '-r300');

%% ================== 3. 图 (b)：目标区放大 ==================
figure('Color','w','Position',[60 60 600 620]);
hold on; grid on; axis equal; view([-42, 24]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

% 真目标圆柱（半径 7 m、高 10 m，下底面圆心 (0,200,0)）
theta = linspace(0, 2*pi, 80);
x_cyl = true_c(1) + Rc*cos(theta);
y_cyl = true_c(2) + Rc*sin(theta);
z_bot = true_c(3);  z_top = true_c(3) + Hc;
% 侧面
Xc = repmat(x_cyl, 2, 1);
Yc = repmat(y_cyl, 2, 1);
Zc = [z_bot*ones(size(x_cyl)); z_top*ones(size(x_cyl))];
surf(Xc, Yc, Zc, 'FaceColor', c_true, 'FaceAlpha', 0.45, 'EdgeColor', 'none');
% 上下底面
fill3(x_cyl, y_cyl, z_top*ones(size(x_cyl)), c_true, 'FaceAlpha', 0.45, 'EdgeColor', 'none');
fill3(x_cyl, y_cyl, z_bot*ones(size(x_cyl)), c_true, 'FaceAlpha', 0.45, 'EdgeColor', 'none');
% 真目标标签：往 +X、+Z 方向移开一点
tc_lab_off = [10; 0; 5];
text(true_c(1)+tc_lab_off(1), true_c(2)+tc_lab_off(2), z_top+3+tc_lab_off(3), '真目标', ...
     'Color', c_true, 'FontSize', fs);

% 假目标（原点）：蓝色圆点
scatter3(0, 0, 0, 100, c_fake, 'o', 'filled');
text(0, -18, 0, '假目标', 'Color', c_fake, 'FontSize', fs);

% 假目标与真目标圆心连线（体现 200 m 偏移）
plot3([0 true_c(1)], [0 true_c(2)], [0 true_c(3)], 'k--', 'LineWidth', 1.0);
text(12, 100, 2, '200 m', 'Color', [0.3 0.3 0.3], 'FontSize', fs-1, ...
     'HorizontalAlignment', 'center');

% 来袭方向示意（导弹从 +X 高空沿斜线直指假目标原点）
plot3([35 3], [0 0], [12 0], '--', 'Color', c_m, 'LineWidth', 1.4);
text(35, 0, 15, '来袭方向', 'Color', c_m, 'FontSize', fs-1);

set(gca, 'FontName', 'Times New Roman', 'FontSize', fs, 'LineWidth', 0.8);
xlim([-40 40]); ylim([-40 240]); zlim([-5 30]);
box on;

%% ================== 4. 输出 ==================
print(gcf, 'background_scene_target', '-dpng', '-r300');
fprintf('已输出：background_scene_overview.png（全局态势）、background_scene_target.png（目标区放大）\n');
