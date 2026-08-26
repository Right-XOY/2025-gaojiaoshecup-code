clear;
clc;
close all;

%% ==================== 场景常量 ====================
params.g = 9.8;
params.R = 10;
params.v_m = 300;
params.M0 = [20000; 0; 2000];
params.u_m = -params.M0 / norm(params.M0);
params.T_end = norm(params.M0) / params.v_m;

params.uav0 = [17800; 0; 1800];
params.cylCenter = [0, 200];
params.cylR = 7;
params.cylTop = 10;
params.cylBottom = 0;

%% ==================== Q2 最优结果 ====================
theta_deg = 8.15;
v_opt = 86.40;
tl_opt = 0.0010;
tau_opt = 0.9760;
td_opt = tl_opt + tau_opt;

% 实际投放点、起爆点：始终用于真实轨迹计算
Pl = [17800.09; 0.01; 1800.00];
Pd = [17883.56; 11.96; 1795.33];

u_h = [cosd(theta_deg); sind(theta_deg); 0];

%% ==================== 仅主图使用的合并显示点 ====================
Pld_mid = (Pl + Pd) / 2;
targetPoint = [0; 0; 0];

%% ==================== 轨迹计算 ====================
nBomb = 150;
t_bomb = linspace(0, tau_opt, nBomb);
bomb_traj = zeros(3, nBomb);

for k = 1:nBomb
    dtb = t_bomb(k);
    bomb_traj(:, k) = Pl ...
        + v_opt * dtb * u_h ...
        - [0; 0; 0.5 * params.g * dtb^2];
end

nCloud = 120;
cdt = linspace(td_opt, min(td_opt + 20, params.T_end), nCloud);
Ctraj = Pd - [0; 0; 3] * (cdt - td_opt);

nMissile = 400;
mt = linspace(0, params.T_end, nMissile);
Mtraj = params.M0' + params.v_m * mt' * params.u_m';

nUAV = 120;
ut = linspace(0, tl_opt, nUAV);
UAVtraj = params.uav0' + v_opt * ut' * u_h';

%% ==================== 真目标圆柱 ====================
theta_c = linspace(0, 2 * pi, 80);
x_cyl = params.cylCenter(1) + params.cylR * cos(theta_c);
y_cyl = params.cylCenter(2) + params.cylR * sin(theta_c);

Xc = repmat(x_cyl, 2, 1);
Yc = repmat(y_cyl, 2, 1);
Zc = [params.cylBottom * ones(size(x_cyl)); ...
      params.cylTop * ones(size(x_cyl))];

%% ==================== 颜色 ====================
colM = [0.82, 0.10, 0.12];
colU = [0.05, 0.25, 0.72];
colB = [0.00, 0.58, 0.68];
colC = [0.72, 0.08, 0.58];
colP = [0.08, 0.58, 0.28];
colTarget = [0.95, 0.42, 0.05];
colDark = [0.12, 0.12, 0.12];
colGold = [1.00, 0.72, 0.05];

%% ==================== 主图 ====================
fig = figure('Color', 'w', 'Position', [100, 100, 900, 680]);
ax_main = axes('Parent', fig);
hold(ax_main, 'on');

axis(ax_main, 'equal');
xlim(ax_main, [-1000, 21000]);
ylim(ax_main, [-1000, 2500]);
zlim(ax_main, [0, 2200]);

set(ax_main, 'FontName', 'Times New Roman', ...
    'FontSize', 10.5, 'LineWidth', 1.0, ...
    'Projection', 'perspective');

xlabel(ax_main, '$X$ (m)', 'Interpreter', 'latex', 'FontSize', 10.5);
ylabel(ax_main, '$Y$ (m)', 'Interpreter', 'latex', 'FontSize', 10.5);
zlabel(ax_main, '$Z$ (m)', 'Interpreter', 'latex', 'FontSize', 10.5);
box(ax_main, 'on');
grid(ax_main, 'off');

% 轨迹
hM = plot3(ax_main, Mtraj(:, 1), Mtraj(:, 2), Mtraj(:, 3), '-', ...
    'Color', colM, 'LineWidth', 1.8);
hU = plot3(ax_main, UAVtraj(:, 1), UAVtraj(:, 2), UAVtraj(:, 3), '-', ...
    'Color', colU, 'LineWidth', 1.8);
hB = plot3(ax_main, bomb_traj(1, :), bomb_traj(2, :), bomb_traj(3, :), '-.', ...
    'Color', colB, 'LineWidth', 1.6);
hC = plot3(ax_main, Ctraj(1, :), Ctraj(2, :), Ctraj(3, :), '-', ...
    'Color', colC, 'LineWidth', 1.8);

% 真目标圆柱仅作为几何轮廓保留
surf(ax_main, Xc, Yc, Zc, ...
    'FaceColor', [0.18, 0.68, 0.30], ...
    'FaceAlpha', 0.28, ...
    'EdgeColor', [0.04, 0.38, 0.14], ...
    'LineWidth', 0.8, ...
    'HandleVisibility', 'off');

% 局部区域框
bbox_x = [17760, 17910, 17910, 17760, 17760];
bbox_y = [-3, -3, 18, 18, -3];
bbox_z = [1772, 1772, 1808, 1808, 1772];
plot3(ax_main, bbox_x, bbox_y, bbox_z, 'r--', ...
    'LineWidth', 1.2, 'HandleVisibility', 'off');

% 起点
sM = scatter3(ax_main, params.M0(1), params.M0(2), params.M0(3), ...
    180, '^', 'filled', 'MarkerFaceColor', colM, ...
    'MarkerEdgeColor', 'w', 'LineWidth', 1.2);
sU = scatter3(ax_main, params.uav0(1), params.uav0(2), params.uav0(3), ...
    150, 'o', 'filled', 'MarkerFaceColor', colU, ...
    'MarkerEdgeColor', 'w', 'LineWidth', 1.2);

% 主图中投放点与起爆点合并显示
hPld = scatter3(ax_main, Pld_mid(1), Pld_mid(2), Pld_mid(3), ...
    250, 'p', 'filled', 'MarkerFaceColor', colP, ...
    'MarkerEdgeColor', colGold, 'LineWidth', 1.6);

% 主图中真、假目标合并为目标点
hTargetPoint = scatter3(ax_main, targetPoint(1), targetPoint(2), targetPoint(3), ...
    260, 'o', 'filled', 'MarkerFaceColor', colTarget, ...
    'MarkerEdgeColor', 'w', 'LineWidth', 1.5);

% 主图文字：白色底色避免和轨迹混叠
text(ax_main, params.M0(1) - 1800, params.M0(2) + 260, params.M0(3) + 80, ...
    'M1', 'Color', colM, 'FontSize', 12, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'Margin', 3, ...
    'EdgeColor', [0.85, 0.55, 0.55], 'Clipping', 'off');

text(ax_main, params.uav0(1) - 1450, params.uav0(2) - 420, params.uav0(3) + 100, ...
    'FY1', 'Color', colU, 'FontSize', 11, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'Margin', 3, ...
    'EdgeColor', [0.55, 0.65, 0.90], 'Clipping', 'off');

text(ax_main, Pld_mid(1) - 2400, Pld_mid(2) + 320, Pld_mid(3) - 180, ...
    '投放点/起爆点', 'Color', [0.03, 0.38, 0.16], ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'Margin', 3, ...
    'EdgeColor', [0.55, 0.80, 0.60], 'Clipping', 'off');

text(ax_main, targetPoint(1) + 850, targetPoint(2) + 520, targetPoint(3) + 300, ...
    '目标点 (0,0,0)', 'Color', [0.85, 0.25, 0.02], ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'Margin', 3, ...
    'EdgeColor', [0.95, 0.65, 0.35], 'Clipping', 'off');

legend(ax_main, ...
    [hM, sM, hU, sU, hB, hPld, hC, hTargetPoint], ...
    {'导弹轨迹', 'M1', '无人机轨迹', 'FY1', ...
     '干扰弹抛射弹道', '投放点/起爆点', ...
     '云团下沉轨迹', '目标点'}, ...
    'Location', 'southwest', 'Box', 'off', 'FontSize', 9.5);

view(ax_main, 3);

%% ==================== 右上角局部放大图 ====================
% 局部图保留真实的投放点 Pl 和起爆点 Pd，不进行合并
ax_inset = axes('Parent', fig, 'Position', [0.61, 0.64, 0.33, 0.26]);
hold(ax_inset, 'on');
axis(ax_inset, 'equal');

xlim(ax_inset, [17750, 17920]);
ylim(ax_inset, [-5, 20]);
zlim(ax_inset, [1770, 1810]);

set(ax_inset, 'FontName', 'Times New Roman', ...
    'FontSize', 8.2, 'LineWidth', 0.8);
title(ax_inset, '投放/起爆局部放大', 'FontSize', 9, 'FontWeight', 'bold');
xlabel(ax_inset, '$X$ (m)', 'Interpreter', 'latex', 'FontSize', 8.2);
ylabel(ax_inset, '$Y$ (m)', 'Interpreter', 'latex', 'FontSize', 8.2);
zlabel(ax_inset, '$Z$ (m)', 'Interpreter', 'latex', 'FontSize', 8.2);
box(ax_inset, 'on');
grid(ax_inset, 'off');

% 局部图轨迹
plot3(ax_inset, bomb_traj(1, :), bomb_traj(2, :), bomb_traj(3, :), '-.', ...
    'Color', colB, 'LineWidth', 1.5);
plot3(ax_inset, Ctraj(1, :), Ctraj(2, :), Ctraj(3, :), '-', ...
    'Color', colC, 'LineWidth', 1.6);

% 局部图保留两个实际点
scatter3(ax_inset, Pl(1), Pl(2), Pl(3), ...
    140, 'd', 'filled', ...
    'MarkerFaceColor', colP, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.0);

scatter3(ax_inset, Pd(1), Pd(2), Pd(3), ...
    160, 'p', 'filled', ...
    'MarkerFaceColor', colDark, ...
    'MarkerEdgeColor', colGold, ...
    'LineWidth', 1.2);

% 局部图标签分别错开，避免相互遮挡
text(ax_inset, Pl(1) - 34, Pl(2) - 2.5, Pl(3) + 5.5, ...
    '投放点', 'Color', [0.03, 0.38, 0.16], ...
    'FontSize', 8, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'Margin', 1, ...
    'EdgeColor', [0.55, 0.80, 0.60], 'Clipping', 'off');

text(ax_inset, Pd(1) + 7, Pd(2) + 1.6, Pd(3) - 3.5, ...
    '起爆点', 'Color', colDark, ...
    'FontSize', 8, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'Margin', 1, ...
    'EdgeColor', [0.90, 0.72, 0.25], 'Clipping', 'off');

plot3(ax_inset, bbox_x, bbox_y, bbox_z, 'r--', ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');

view(ax_inset, 3);

%% ==================== 导出图片 ====================
try
    exportgraphics(fig, 'q2_3d_inset.png', ...
        'Resolution', 300);
catch
    print(fig, 'q2_3d_inset', '-dpng', '-r300');
end

fprintf('绘图完成：q2_3d_inset.png\n');