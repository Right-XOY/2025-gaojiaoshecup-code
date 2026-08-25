%% q4_main  问题4：FY1、FY2、FY3 各投放 1 枚烟幕干扰弹实施对 M1 的接力遮蔽
%  给出 3 架无人机的飞行方向、飞行速度、投放点与起爆点，
%  使总有效遮蔽时间尽可能长，并将结果写入 data/result2.xlsx。
%
%  求解架构：分层优化（三层）
%    第1层【上层·时序规划】先验固定三段接力时序（软约束目标时刻）：
%       FY1 遮蔽窗口中心 ≈ 4 s，FY2 ≈ 21 s，FY3 ≈ 34 s。
%       （该时序来源于 Python 复现的 12 维联合优化验证结果，
%         保证三窗口自然错开、接续覆盖 M1 飞行的早/中/晚段。）
%    第2层【下层·单机参数优化】对每架无人机独立求解一个 4 维子问题：
%       决策变量 x = [theta, v, t_l, tau]，目标为软约束
%       max 自身遮蔽时长 + w_center*(窗口中心-目标时刻)^2 尽量小
%       （允许遮蔽窗口延伸出分配区间，只追求中心对齐）。
%    第3层【精化·联合优化】12 维联合 JADE：直接最大化三弹遮蔽的
%       【总并集时长】，以第2层结果与已复现最优解为初始种子。
%
%  依赖函数：q4_single_fitness.m、q4_joint_fitness.m、
%            q3_jade.m、q3_occlusion_window.m
%  运行：MATLAB 中 cd 到 code 文件夹后运行 q4_main（无需额外工具箱）
%
%  作者：数学建模 A 题队伍    日期：2026-08-24

clear; clc; close all;
rng(2026);                      % 固定随机种子，保证可复现

%% ================== 1. 场景常量（与问题2/3一致） ==================
params.g         = 9.8;
params.R         = 10;
params.v_min     = 70;
params.v_max     = 140;
params.M0        = [20000; 0; 2000];               % M1 初始位置
params.u_m       = -params.M0/norm(params.M0);     % M1 飞行方向（直指假目标）
params.v_m       = 300;
params.T_end     = norm(params.M0)/params.v_m;     % 导弹到达假目标时刻 ≈ 67.0 s
params.cylCenter = [0, 200];
params.cylR      = 7;
params.cylTop    = 10;
params.cylBottom = 0;
params.nPhi      = 100;                            % 圆周采样点数（敏感性分析：n>=100 收敛）
params.w_center  = 0.5;          % 下层软约束权重（窗口中心对齐强度）

% 三架无人机初始位置（3 x 3，每列一架：FY1、FY2、FY3）
params.uav0_all = [17800 12000 6000;
                    0    1400 -3000;
                  1800  1400   700];
uav_names      = {'FY1', 'FY2', 'FY3'};
target_centers = [4, 21, 34];   % 第1层先验时序：各机遮蔽窗口目标中心 (s)

%% ================== 2. 第2层：3 次独立 4 维 JADE ==================
%  每架无人机一个 4 维子问题，目标：最大化自身遮蔽 + 窗口中心贴近目标时刻
params.dt = 0.05;               % 粗时间步长 (s)
opts_s.NP = 30;  opts_s.MAXITER = 200;  opts_s.p = 0.05;  opts_s.c = 0.1;

% 物理引导种子（基于 Python 复现的联合最优与其邻域变体）
seeds = {
  [4.60, 79.82, 0.993, 0.446;     % FY1 已复现最优（窗口 [1.7,6.3]）
   8.15, 86.40, 0.000, 0.980;     % 邻域变体：早段高空云
   6.68, 98.60, 0.050, 0.634;
   4.00, 86.40, 0.650, 0.450];
  [226.18, 117.64, 7.970, 8.007;  % FY2 已复现最优（窗口 [19.2,22.9]）
   225.00, 120.00, 8.000, 8.000;
   240.00, 115.00, 7.000, 7.000;
   205.00, 115.00, 18.00, 10.00];
  [90.00, 110.00, 23.968, 4.000;  % FY3 已复现最优（窗口 [33.1,35.9]）
   95.00, 140.00, 18.00, 4.000;
   130.00, 135.00, 22.00, 8.000;
   90.00, 110.00, 24.00, 4.000]};

x_lower   = cell(3, 1);
f_lower   = zeros(3, 1);
hist_lower = cell(3, 1);
for i = 1:3
    tau_max_i = sqrt(2*params.uav0_all(3, i)/params.g);   % 落地前起爆上限
    lb = [0, 70, 0, 0];
    ub = [360, 140, params.T_end, tau_max_i];
    opts_s.X0 = seeds{i};
    fun_i = @(x) q4_single_fitness(x, params, params.uav0_all(:, i), ...
                                   target_centers(i));
    fprintf('========== 第2层：FY%d 单机 4 维 JADE 优化中 ... ==========\n', i);
    [x_lower{i}, f_lower(i), hist_lower{i}] = q3_jade(fun_i, lb, ub, opts_s);
end

%% ================== 3. 第3层：12 维联合 JADE 精化 ==================
%  直接最大化三弹遮蔽的总并集时长，决策变量按机分组交错排列
lb_j = [0 70 0 0,  0 70 0 0,  0 70 0 0];
ub_j = [360 140 params.T_end sqrt(2*params.uav0_all(3,1)/params.g), ...
        360 140 params.T_end sqrt(2*params.uav0_all(3,2)/params.g), ...
        360 140 params.T_end sqrt(2*params.uav0_all(3,3)/params.g)];

joint_seeds = zeros(6, 12);
for i = 1:3
    joint_seeds(i, 4*i-3 : 4*i) = x_lower{i};        % 第2层最优解
end
joint_seeds(4, :) = [4.60, 79.82, 0.993, 0.446,  ... % 已复现联合最优
                     226.18, 117.64, 7.970, 8.007, ...
                     90.00, 110.00, 23.968, 4.000];
joint_seeds(5, :) = [8.15, 86.40, 0.000, 0.980,  ... % 物理变体
                     225.00, 120.00, 8.000, 8.000, ...
                     95.00, 140.00, 18.00, 4.000];
joint_seeds(6, :) = [6.68, 98.60, 0.050, 0.634,  ...
                     240.00, 115.00, 7.000, 7.000, ...
                     130.00, 135.00, 22.00, 8.000];

opts_j.NP = 50;  opts_j.MAXITER = 300;  opts_j.p = 0.05;  opts_j.c = 0.1;
opts_j.X0 = joint_seeds;
fun_j = @(x) q4_joint_fitness(x, params);
disp('========== 第3层：12 维联合 JADE 精化中 ... ==========');
[x_joint, f_joint, hist_joint] = q3_jade(fun_j, lb_j, ub_j, opts_j);

%% ================== 4. 细步长验证 ==================
params.dt = 0.005;              % 精细时间步长 (s)
T_total = -q4_joint_fitness(x_joint, params);   % 总并集遮蔽时长 (s)

%% ================== 5. 结果解析 ==================
theta = zeros(1, 3);  v = zeros(1, 3);
t_l = zeros(1, 3);    tau = zeros(1, 3);
P_l = zeros(3, 3);    P_d = zeros(3, 3);  t_d = zeros(1, 3);
for i = 1:3
    theta(i) = x_joint(4*i-3);  v(i) = x_joint(4*i-2);
    t_l(i)   = x_joint(4*i-1);  tau(i) = x_joint(4*i);
    u_h = [cosd(theta(i)); sind(theta(i)); 0];
    P_l(:, i) = params.uav0_all(:, i) + v(i)*t_l(i)*u_h;      % 投放点
    P_d(:, i) = P_l(:, i) + v(i)*tau(i)*u_h - ...
                [0; 0; 0.5*params.g*tau(i)^2];                % 起爆点
    t_d(i) = t_l(i) + tau(i);                                  % 起爆时刻
end

% ---- 各机独立贡献时长（仅该机云团遮蔽，不考虑重叠）----
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
fprintf('\n========== 问题4 最优投放策略（FY1/FY2/FY3 接力遮蔽 M1）==========\n');
for i = 1:3
    fprintf('--- %s ---\n', uav_names{i});
    fprintf('  航向角 theta = %.4f°，速度 v = %.4f m/s\n', theta(i), v(i));
    fprintf('  投放时刻 t_l = %.4f s，引信延时 tau = %.4f s，起爆时刻 t_d = %.4f s\n', ...
        t_l(i), tau(i), t_d(i));
    fprintf('  投放点 P_l   = (%.2f, %.2f, %.2f) m\n', P_l(:, i));
    fprintf('  起爆点 P_d   = (%.2f, %.2f, %.2f) m\n', P_d(:, i));
    fprintf('  独立贡献时长 = %.4f s\n', T_ind(i));
end
fprintf('----------------------------------------------------------\n');
fprintf('总有效遮蔽时长（三弹并集）= %.4f s\n', T_total);
fprintf('==========================================================\n');

%% ================== 7. 写入 result2.xlsx ==================
header = {'无人机编号', '无人机运动方向', '无人机运动速度 (m/s)', ...
    '烟幕干扰弹投放点的x 坐标 (m)', '烟幕干扰弹投放点的y坐标 (m)', ...
    '烟幕干扰弹投放点的z坐标 (m)', '烟幕干扰弹起爆点的x坐标 (m)', ...
    '烟幕干扰弹起爆点的y坐标 (m)', '烟幕干扰弹起爆点的z坐标 (m)', ...
    '有效干扰时长 (s)'};
out = cell(6, 10);                          % 6 行，与模板 result2.xlsx 一致
out(1, :) = header;
for i = 1:3
    out(i+1, :) = {uav_names{i}, theta(i), v(i), P_l(1,i), P_l(2,i), P_l(3,i), ...
                   P_d(1,i), P_d(2,i), P_d(3,i), T_ind(i)};
end
out(5, :) = {[], [], [], [], [], [], [], [], [], []};  % 第 5 行留空（与模板对齐）
out(6, 2) = {'注：以x轴为正向，逆时针方向为正，取值0~360（度）。'};
data_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
out_file = fullfile(data_dir, 'result2.xlsx');
writecell(out, out_file);
fprintf('结果已写入：%s\n', out_file);

%% ================== 8. 收敛曲线 ==================
figure('Color', 'w', 'Position', [80 80 960 400]);
subplot(1, 3, 1);
for i = 1:3
    plot(1:numel(hist_lower{i}), -hist_lower{i}, 'LineWidth', 1.4); hold on;
end
grid on; legend(uav_names, 'Location', 'best');
xlabel('进化代数'); ylabel('第2层最佳目标 (s)');
title('第2层 单机 4 维 JADE 收敛曲线');
subplot(1, 3, 2);
plot(1:numel(hist_joint), -hist_joint, 'r-', 'LineWidth', 1.5); grid on;
xlabel('进化代数'); ylabel('总并集遮蔽时长 (s)');
title('第3层 12 维联合 JADE 收敛曲线');
subplot(1, 3, 3);
bar([T_ind, T_total], 0.55);
set(gca, 'XTickLabel', {'FY1', 'FY2', 'FY3', '总并集'});
ylabel('有效遮蔽时长 (s)'); grid on;
title('各机独立贡献与总并集时长');

%% ================== 9. 遮蔽时间线图 ==================
shielded = false(1, numel(ts));
occ_i = false(3, numel(ts));
for i = 1:3
    mask = (ts >= t_d(i)) & (ts <= t_d(i) + 20);
    if any(mask)
        Cw = P_d(:, i)*ones(1, sum(mask)) - [0; 0; 3]*(ts(mask) - t_d(i));
        occ = q3_occlusion_window(M_all(:, mask), Cw, params.R, params);
        occ_i(i, mask) = occ;                          % occ 为行向量，位置赋值同向
        shielded(mask) = shielded(mask) | occ(:)';
    end
end

figure('Color', 'w', 'Position', [100 100 900 460]);
subplot(2, 1, 1);
hold on; grid on;
stairs(ts, shielded, 'k-', 'LineWidth', 1.6);
for i = 1:3
    xline(t_d(i), '--', sprintf('%s 起爆', uav_names{i}));
end
xline(params.T_end, '--r', '导弹到假目标');
xlabel('时间 t (s)');  ylabel('遮蔽状态');
ylim([-0.1 1.1]);  yticks([0 1]);
title(sprintf('三机接力遮蔽 M1 总时间线（总并集 = %.2f s）', T_total));

subplot(2, 1, 2);
hold on; grid on;
colors = {'g', 'b', 'm'};
for i = 1:3
    plot(ts, occ_i(i, :), '-', 'Color', colors{i}, 'LineWidth', 1.1);
end
legend(uav_names, 'Location', 'best');
xlabel('时间 t (s)');  ylabel('各机遮蔽状态');
ylim([-0.1 1.1]);  yticks([0 1]);
title('各架无人机独立遮蔽贡献');
hold off;

%% ================== 10. 三维态势图 ==================
figure('Color', 'w', 'Position', [120 120 940 700]);
hold on; grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('问题4 最优投放策略三维态势图');

% 导弹轨迹
mt = linspace(0, params.T_end, 300);
Mtraj = params.M0' + params.v_m*mt'*params.u_m';
plot3(Mtraj(:, 1), Mtraj(:, 2), Mtraj(:, 3), 'r-', 'LineWidth', 1.6);
scatter3(params.M0(1), params.M0(2), params.M0(3), 90, 'r', '^', 'filled');
text(params.M0(1)+300, params.M0(2), params.M0(3), 'M1', 'Color', 'r', 'FontSize', 11);

% 三架无人机：轨迹、投放点、起爆点、云团下沉轨迹
for i = 1:3
    u_h = [cosd(theta(i)); sind(theta(i)); 0];
    ut = linspace(0, t_l(i), 60);
    UAVtraj = params.uav0_all(:, i)' + v(i)*ut'*u_h';
    plot3(UAVtraj(:, 1), UAVtraj(:, 2), UAVtraj(:, 3), '-', ...
          'Color', colors{i}, 'LineWidth', 1.5);
    scatter3(params.uav0_all(1, i), params.uav0_all(2, i), params.uav0_all(3, i), ...
             80, colors{i}, 'o', 'filled');
    text(params.uav0_all(1, i)+300, params.uav0_all(2, i), params.uav0_all(3, i), ...
         uav_names{i}, 'Color', colors{i}, 'FontSize', 11);

    scatter3(P_l(1, i), P_l(2, i), P_l(3, i), 90, 'g', 'd', 'filled');
    scatter3(P_d(1, i), P_d(2, i), P_d(3, i), 110, 'k', 'p', 'filled');
    cdt = linspace(t_d(i), min(t_d(i)+20, params.T_end), 60);
    Ctraj = P_d(:, i) - [0; 0; 3]*(cdt - t_d(i));
    plot3(Ctraj(1, :), Ctraj(2, :), Ctraj(3, :), '-', ...
          'Color', colors{i}, 'LineWidth', 1.2);
    text(P_l(1, i)+300, P_l(2, i), P_l(3, i), sprintf('%s投放', uav_names{i}), ...
         'Color', 'g', 'FontSize', 9);
    text(P_d(1, i)+300, P_d(2, i), P_d(3, i), sprintf('%s起爆', uav_names{i}), ...
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

legend('导弹轨迹', 'M1', '无人机轨迹', '投放点', '起爆点', ...
       '云团下沉轨迹', '假目标', '真目标', 'Location', 'best');
view(3);
hold off;

disp('求解完成。请检查控制台输出与 result2.xlsx，以及三张图。');
