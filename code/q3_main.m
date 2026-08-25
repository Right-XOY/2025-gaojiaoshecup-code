%% q3_main  问题3：利用无人机 FY1 投放 3 枚烟幕干扰弹实施对 M1 的干扰
%  给出 FY1 的飞行方向、飞行速度、3 枚烟幕干扰弹的投放点与起爆点，
%  使总有效遮蔽时间尽可能长，并将结果写入 data/result1.xlsx。
%
%  求解算法：JADE 自适应差分进化（带外部存档）
%    - 全局搜索：粗时间步长 JADE（dt=0.1 s）
%    - 精化验证：细时间步长 JADE（dt=0.005 s）以粗解为初始点继续进化
%    - 约束（投放间隔>=1s、起爆在交战窗口内）以惩罚函数处理
%
%  依赖函数：q3_fitness.m、q3_jade.m、q3_occlusion_window.m
%  运行：MATLAB 中 cd 到 code 文件夹后运行 q3_main（无需额外工具箱）
%
%  作者：数学建模 A 题队伍    日期：2026-08-24

clear; clc; close all;
rng(2026);                      % 固定随机种子，保证可复现

%% ================== 1. 场景常量（与问题2一致） ==================
params.g         = 9.8;
params.R         = 10;
params.v_min     = 70;
params.v_max     = 140;
params.uav0      = [17800; 0; 1800];               % FY1 初始位置
params.M0        = [20000; 0; 2000];               % M1 初始位置
params.u_m       = -params.M0/norm(params.M0);     % M1 飞行方向（直指假目标）
params.v_m       = 300;
params.T_end     = norm(params.M0)/params.v_m;     % 导弹到达假目标时刻 ≈ 67.0 s
params.cylCenter = [0, 200];
params.cylR      = 7;
params.cylTop    = 10;
params.cylBottom = 0;
params.nPhi      = 100;                            % 圆周采样点数（敏感性分析：n>=100 收敛）
params.tau_max   = sqrt(2*1800/params.g);          % 引信延时上限 ≈ 19.2 s

%% ================== 2. 决策变量与边界 ==================
% x = [theta_deg, v, t_l1, tau1, t_l2, tau2, t_l3, tau3]
nvars = 8;
lb = [0,    70,  0, 0,  0, 0,  0, 0];
ub = [360,  140, params.T_end, params.tau_max, ...
              params.T_end, params.tau_max, params.T_end, params.tau_max];

%% ================== 3. 粗步长 JADE 全局搜索 ==================
params.dt = 0.1;                % 粗时间步长 (s)
opts.NP = 60;  opts.MAXITER = 300;  opts.p = 0.05;  opts.c = 0.1;

% 基于物理预扫描的初始种群（帮助 JADE 进入有利盆域）：
%   x = [theta_deg, v, t_l1, tau1, t_l2, tau2, t_l3, tau3]
%   盆域A：-x 链式（θ≈180°，三弹沿 -x 依次投放，窗口 [5,12.6] 连续衔接，
%          经 Python 复现预扫描确认的全局最优区，总遮蔽约 7.4 s）
%   盆域B：+x 早段（θ≈0~10°，云团在 FY1 起始区高空，窗口 [1.5,7.3]）
%   其余：覆盖全空间随机散布
init_pop = [
    179.64, 137.74, 0.008, 3.672, 3.287, 5.043, 5.383, 5.882;  % A 已知最优附近
    180.00, 140.00, 0.000, 3.500, 3.100, 5.300, 5.200, 6.100;
    179.50, 135.00, 0.400, 3.800, 1.400, 4.300, 4.300, 5.500;
    179.50, 130.00, 0.400, 3.700, 3.900, 5.200, 5.200, 5.700;
    179.50, 140.00, 0.200, 3.700, 1.200, 4.300, 4.000, 5.500;
    180.00, 135.00, 0.000, 3.400, 3.100, 5.200, 5.200, 5.900;
    180.00, 130.00, 0.000, 3.300, 3.100, 5.100, 5.200, 5.700;
    179.00, 125.00, 0.200, 3.400, 1.200, 4.200, 2.600, 5.000;
    4.00,   86.40, 0.650, 0.450, 1.750, 0.050, 2.750, 0.050;   % B 早段
    6.68,   98.60, 0.050, 0.634, 1.100, 0.407, 2.500, 0.500;
    8.15,   86.40, 0.000, 0.980, 1.500, 1.000, 2.500, 1.000;
    0.00,  100.00, 1.000, 0.600, 2.000, 0.400, 3.000, 0.500;   % 随机散布
    90.00, 120.00, 2.000, 3.000, 5.000, 4.000, 8.000, 5.000;
    270.00,110.00, 1.000, 4.000, 4.000, 5.000, 8.000, 6.000;
    45.00,  90.00, 5.000, 3.000, 9.000, 4.000, 14.00, 5.000;
    300.00,130.00, 3.000, 4.500, 7.000, 5.500, 12.00, 6.500;
    150.00,110.00, 6.000, 4.000, 10.00, 5.000, 15.00, 6.000;
    20.00, 120.00, 0.500, 2.000, 2.500, 3.000, 4.500, 4.000];
opts.X0 = init_pop;

fun = @(x) q3_fitness(x, params);
disp('========== JADE 粗步长全局搜索中 ... ==========');
[x_opt, f_opt, hist_best] = q3_jade(fun, lb, ub, opts);
disp(sprintf('粗解总遮蔽时长 = %.4f s', -f_opt));

%% ================== 4. 细步长 JADE 精化 ==================
params.dt = 0.005;              % 精细时间步长 (s)
opts_fine.NP = 20;  opts_fine.MAXITER = 120;
opts_fine.p = 0.05;  opts_fine.c = 0.1;
opts_fine.X0 = x_opt;           % 以粗解作为初始种群种子
fun_fine = @(x) q3_fitness(x, params);
disp('========== JADE 细步长精化中 ... ==========');
[x_opt, f_opt, hist_best_fine] = q3_jade(fun_fine, lb, ub, opts_fine);
T_total = -f_opt;               % 总并集遮蔽时长 (s)

%% ================== 5. 结果解析 ==================
theta_deg = x_opt(1);  v = x_opt(2);
t_l = x_opt(3:2:7);    tau = x_opt(4:2:8);
u_h = [cosd(theta_deg); sind(theta_deg); 0];
P_l = params.uav0 + v*t_l.*u_h;                    % 3x3：各弹投放点
P_d = P_l + v*tau.*u_h - [zeros(1,3); zeros(1,3); 0.5*params.g*tau.^2];  % 3x3：各弹起爆点
t_d = t_l + tau;                                   % 各弹起爆时刻

% ---- 各弹独立贡献时长（仅该弹云团遮蔽的时间，不考虑重叠）----
ts = 0 : params.dt : params.T_end;
M_all = params.M0*ones(1, numel(ts)) + params.v_m*(params.u_m*ts);
T_ind = zeros(1, 3);
for i = 1:3
    idx = find(ts >= t_d(i) & ts <= t_d(i) + 20);
    if ~isempty(idx)
        Cw = P_d(:, i)*ones(1, numel(idx)) - [0; 0; 3]*(ts(idx) - t_d(i));
        T_ind(i) = sum(q3_occlusion_window(M_all(:, idx), Cw, params.R, params))*params.dt;
    end
end

%% ================== 6. 控制台输出 ==================
fprintf('\n========== 问题3 最优投放策略（FY1 投放 3 枚干扰弹干扰 M1）==========\n');
fprintf('飞行航向角 theta = %.4f°\n', theta_deg);
fprintf('飞行速度   v     = %.4f m/s\n', v);
for i = 1:3
    fprintf('--- 第 %d 枚弹 ---\n', i);
    fprintf('  投放时刻 t_l = %.4f s，引信延时 tau = %.4f s，起爆时刻 t_d = %.4f s\n', ...
        t_l(i), tau(i), t_d(i));
    fprintf('  投放点 P_l   = (%.2f, %.2f, %.2f) m\n', P_l(:, i));
    fprintf('  起爆点 P_d   = (%.2f, %.2f, %.2f) m\n', P_d(:, i));
    fprintf('  独立贡献时长 = %.4f s\n', T_ind(i));
end
fprintf('----------------------------------------------------------\n');
fprintf('总有效遮蔽时长（三弹并集）= %.4f s\n', T_total);
fprintf('==========================================================\n');

%% ================== 7. 写入 result1.xlsx ==================
header = {'无人机运动方向', '无人机运动速度 (m/s)', '烟幕干扰弹编号', ...
    '烟幕干扰弹投放点的x坐标 (m)', '烟幕干扰弹投放点的y坐标 (m)', ...
    '烟幕干扰弹投放点的z坐标 (m)', '烟幕干扰弹起爆点的x坐标 (m)', ...
    '烟幕干扰弹起爆点的y坐标 (m)', '烟幕干扰弹起爆点的z坐标 (m)', ...
    '有效干扰时长 (s)'};
out = cell(6, 10);                          % 6 行，与模板 result1.xlsx 一致
out(1, :) = header;
for i = 1:3
    out(i+1, :) = {theta_deg, v, i, P_l(1,i), P_l(2,i), P_l(3,i), ...
                   P_d(1,i), P_d(2,i), P_d(3,i), T_ind(i)};
end
out(5, :) = {[], [], [], [], [], [], [], [], [], []};  % 第 5 行留空（与模板对齐）
out(6, 1) = {'注：以x轴为正向，逆时针方向为正，取值0~360（度）。'};
data_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
out_file = fullfile(data_dir, 'result1.xlsx');
writecell(out, out_file);
fprintf('结果已写入：%s\n', out_file);

%% ================== 8. 收敛曲线 ==================
figure('Color', 'w', 'Position', [80 80 900 400]);
subplot(1, 2, 1);
plot(1:numel(hist_best), -hist_best, 'b-', 'LineWidth', 1.5); grid on;
xlabel('进化代数'); ylabel('最佳总遮蔽时长 (s)');
title('粗步长 JADE 收敛曲线');
subplot(1, 2, 2);
plot(1:numel(hist_best_fine), -hist_best_fine, 'r-', 'LineWidth', 1.5); grid on;
xlabel('进化代数'); ylabel('最佳总遮蔽时长 (s)');
title('细步长 JADE 收敛曲线');

%% ================== 9. 遮蔽时间线图 ==================
shielded = false(1, numel(ts));
for i = 1:3
    mask = (ts >= t_d(i)) & (ts <= t_d(i) + 20);
    if any(mask)
        Cw = P_d(:, i)*ones(1, sum(mask)) - [0; 0; 3]*(ts(mask) - t_d(i));
        occ = q3_occlusion_window(M_all(:, mask), Cw, params.R, params);
        shielded(mask) = shielded(mask) | occ(:)';
    end
end

figure('Color', 'w', 'Position', [100 100 860 420]);
hold on; grid on;
edges = diff([0, shielded, 0]);
starts = find(edges == 1)*params.dt - params.dt;
stops  = find(edges == -1)*params.dt - params.dt;
for k = 1:numel(starts)
    fill([starts(k), stops(k), stops(k), starts(k)], ...
         [0, 0, 1, 1], [0.9 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
stairs(ts, shielded, 'k-', 'LineWidth', 1.5);
colors = {'g', 'b', 'm'};
for i = 1:3
    xline(t_d(i), '--', sprintf('弹%d起爆', i), 'Color', colors{i});
end
xline(params.T_end, '--r', '导弹到假目标');
xlabel('时间 t (s)');  ylabel('遮蔽状态');
ylim([-0.15 1.15]);  yticks([0 1]);
title(sprintf('M1 对真目标的总遮蔽时间线（三弹并集总时长 = %.2f s）', T_total));
hold off;

%% ================== 10. 三维态势图 ==================
figure('Color', 'w', 'Position', [120 120 940 700]);
hold on; grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('问题3 最优投放策略三维态势图');

% 导弹轨迹
mt = linspace(0, params.T_end, 300);
Mtraj = params.M0' + params.v_m*mt'*params.u_m';
plot3(Mtraj(:, 1), Mtraj(:, 2), Mtraj(:, 3), 'r-', 'LineWidth', 1.6);
scatter3(params.M0(1), params.M0(2), params.M0(3), 90, 'r', '^', 'filled');
text(params.M0(1)+300, params.M0(2), params.M0(3), 'M1', 'Color', 'r', 'FontSize', 11);

% 无人机轨迹（到最后一枚弹投放为止）
ut = linspace(0, max(t_l), 100);
UAVtraj = params.uav0' + v*ut'*u_h';
plot3(UAVtraj(:, 1), UAVtraj(:, 2), UAVtraj(:, 3), 'b-', 'LineWidth', 1.6);
scatter3(params.uav0(1), params.uav0(2), params.uav0(3), 80, 'b', 'o', 'filled');
text(params.uav0(1)+300, params.uav0(2), params.uav0(3), 'FY1', 'Color', 'b', 'FontSize', 11);

% 三枚弹：投放点、起爆点、云团下沉轨迹
for i = 1:3
    scatter3(P_l(1, i), P_l(2, i), P_l(3, i), 90, 'g', 'd', 'filled');
    scatter3(P_d(1, i), P_d(2, i), P_d(3, i), 110, 'k', 'p', 'filled');
    cdt = linspace(t_d(i), min(t_d(i)+20, params.T_end), 60);
    Ctraj = P_d(:, i) - [0; 0; 3]*(cdt - t_d(i));
    plot3(Ctraj(1, :), Ctraj(2, :), Ctraj(3, :), 'm-', 'LineWidth', 1.4);
    text(P_l(1, i)+300, P_l(2, i), P_l(3, i), sprintf('投放%d', i), ...
         'Color', 'g', 'FontSize', 9);
    text(P_d(1, i)+300, P_d(2, i), P_d(3, i), sprintf('起爆%d', i), ...
         'Color', 'k', 'FontSize', 9);
end

% 假目标与真目标
scatter3(0, 0, 0, 90, 'm', '*');
text(500, 0, -100, '假目标(0,0,0)', 'Color', 'm', 'FontSize', 10);
theta_c = linspace(0, 2*pi, 60);
x_cyl = params.cylCenter(1) + params.cylR*cos(theta_c);
y_cyl = params.cylCenter(2) + params.cylR*sin(theta_c);
surf(repmat(x_cyl, 2, 1), repmat(y_cyl, 2, 1), ...
     [params.cylBottom*ones(size(x_cyl)); params.cylTop*ones(size(x_cyl))], ...
     'FaceColor', 'g', 'FaceAlpha', 0.25, 'EdgeColor', 'g', 'LineWidth', 0.8);
text(params.cylCenter(1), params.cylCenter(2), params.cylTop+15, ...
     '真目标圆柱 r=7m h=10m', 'Color', 'g', 'FontSize', 10);

legend('导弹轨迹', 'M1', '无人机轨迹', 'FY1', '投放点', '起爆点', ...
       '云团下沉轨迹', '假目标', '真目标', 'Location', 'best');
view(3);
hold off;

disp('求解完成。请检查控制台输出与 result1.xlsx，以及三张图。');
