%% q2_main  问题2：利用无人机 FY1 投放 1 枚烟幕干扰弹实施对 M1 的干扰
%  确定 FY1 的飞行方向、飞行速度、烟幕干扰弹投放点、起爆点，使遮蔽时间尽可能长。
%
%  求解算法：文化基因算法（Memetic Algorithm）
%     - 全局搜索：遗传算法 ga（实数编码，自动处理边界与线性约束）
%     - 局部搜索：模式搜索 patternsearch（作为 ga 的 HybridFcn，对精英解精化）
%     - 之后再用细时间步长对最优解做一次 patternsearch 精化与高精度复核
%
%  说明：需要 Global Optimization Toolbox（提供 ga 与 patternsearch）。
%  运行：在 MATLAB 中 cd 到 code 文件夹后直接运行 q2_main。
%
%  作者：数学建模 A 题队伍    日期：2026-08-24

clear; clc; close all;
rng(2025);                      % 固定随机种子，保证结果可复现

%% ================== 1. 场景常量 ==================
params.g         = 9.8;                            % 重力加速度 (m/s^2)
params.R         = 10;                             % 烟幕有效遮蔽半径 (m)
params.v_min     = 70;                             % 无人机速度下限 (m/s)
params.v_max     = 140;                            % 无人机速度上限 (m/s)
params.uav0      = [17800; 0; 1800];               % FY1 初始位置
params.M0        = [20000; 0; 2000];               % M1 初始位置
params.u_m       = -params.M0/norm(params.M0);     % M1 飞行方向（直指假目标原点）
params.v_m       = 300;                            % 导弹速度 (m/s)
params.T_end     = norm(params.M0)/params.v_m;     % 导弹到达假目标时刻 (s) ≈ 67.0
params.cylCenter = [0, 200];                       % 真目标圆柱轴心 (x,y)
params.cylR      = 7;                              % 圆柱半径 (m)
params.cylTop    = 10;                             % 圆柱顶面高度 (m)
params.cylBottom = 0;                              % 圆柱底面高度 (m)
params.nPhi      = 100;                            % 每圆周采样点数（敏感性分析：n>=100 收敛）
params.tau_max   = sqrt(2*1800/params.g);          % 引信延时上限：保证落地前起爆 (s) ≈ 19.2

%% ================== 2. 决策变量与约束 ==================
% x = [theta, v, t_l, tau]
%   theta : 航向角 [0, 2*pi]
%   v     : 速度   [70, 140]
%   t_l   : 投放时刻 [0, T_end]
%   tau   : 引信延时 [0, tau_max]
nvars = 4;
lb = [0,      params.v_min, 0, 0];
ub = [2*pi,   params.v_max, params.T_end, params.tau_max];
Aineq = [0, 0, 1, 1];           % 线性约束：t_l + tau <= T_end（起爆须在导弹到达假目标前）
bineq = params.T_end;

%% ================== 3. 文化基因算法求解 ==================
params.dt = 0.05;               % 遗传/局部搜索阶段粗时间步长 (s)
fun = @(x) q2_fitness(x, params);

% 局部搜索（patternsearch）选项
psopts = optimoptions('patternsearch', 'MeshTolerance', 1e-3, 'Display', 'off');

% 基于物理分析的初始种群（帮助 ga 快速进入有利盆域）：
%   x = [theta, v, t_l, tau]
%   盆域A：立即投放 + 约1s起爆，云团在FY1起始区高空拦截导弹早期视线（经预扫描确认的全局最优区）
%   盆域B：中段投放，云团落在导弹路径 x≈16000~17000 高空
%   其余：覆盖全空间随机散布
init_pop = [
    0.1422,  86.4, 0.001, 0.976;   % A 已知最优附近
    0.1300,  82.0, 0.000, 0.900;
    0.1500,  80.0, 0.000, 0.900;
    0.1200,  85.0, 0.000, 1.000;
    0.1600,  88.0, 0.000, 1.000;
    0.1000,  78.0, 0.500, 0.800;
    0.1400,  84.0, 1.000, 1.100;
    0.1100,  90.0, 0.300, 0.700;
    0.1700,  86.0, 0.200, 1.200;
    pi,      140.0, 4.000, 5.500;   % B 中段投放
    pi,      130.0, 5.000, 6.000;
    pi+0.1,  120.0, 3.000, 5.000;
    pi-0.1,  140.0, 6.000, 6.500;
    3.0000,  120.0, 4.000, 6.000;
    3.3000,  140.0, 5.000, 5.000;
    0.5000,  100.0, 10.0,  5.000;   % 随机散布
    1.0000,  120.0, 15.0,  8.000;
    2.0000,   70.0, 20.0, 10.0;
    2.5000,  110.0,  8.0, 12.0;
    4.0000,   90.0, 25.0,  6.0;
    5.0000,  140.0, 30.0, 10.0;
    0.3000,  130.0, 40.0,  5.0;
    3.5000,  100.0, 50.0, 10.0;
    1.5000,   75.0,  5.0,  3.0;
    2.2000,  135.0, 12.0,  7.0;
    5.5000,   85.0, 20.0,  4.0;
    0.8000,  110.0, 35.0, 10.0];

% 遗传算法选项（含 HybridFcn -> patternsearch，构成文化基因算法）
gaopts = optimoptions('ga', ...
    'PopulationSize',       80, ...
    'InitialPopulationMatrix', init_pop, ...
    'MaxGenerations',       200, ...
    'MaxStallGenerations',  30, ...
    'FunctionTolerance',    1e-4, ...
    'HybridFcn',            {@patternsearch, psopts}, ...
    'Display',              'iter', ...
    'PlotFcn',              {});

disp('========== 文化基因算法求解中（ga + patternsearch）==========');
[x_opt, ~] = ga(fun, nvars, Aineq, bineq, [], [], lb, ub, [], gaopts);

% ---- 细步长精化：多起点局部搜索取最优 ----
params.dt = 0.005;              % 精细时间步长 (s)
fun_fine = @(x) q2_fitness(x, params);
psopts_fine = optimoptions('patternsearch', 'MeshTolerance', 1e-5, 'Display', 'off');

x_candidates = [x_opt; init_pop(1:9, :)];   % ga 最优 + 盆域A关键种子
f_best = inf;
for i = 1:size(x_candidates, 1)
    x_try = patternsearch(fun_fine, x_candidates(i, :), Aineq, bineq, [], [], lb, ub, [], psopts_fine);
    f_try = fun_fine(x_try);
    if f_try < f_best
        f_best = f_try;
        x_opt = x_try;
    end
end
T_shield = -f_best;             % 高精度遮蔽时长 (s)

%% ================== 4. 结果整理 ==================
theta = x_opt(1);  v = x_opt(2);  t_l = x_opt(3);  tau = x_opt(4);
u_h = [cos(theta); sin(theta); 0];
P_l = params.uav0 + v*t_l*u_h;                          % 投放点
P_d = P_l + v*tau*u_h - [0; 0; 0.5*params.g*tau^2];     % 起爆点
t_det = t_l + tau;                                      % 起爆时刻

fprintf('\n========== 问题2 最优投放策略（FY1 干扰 M1）==========\n');
fprintf('飞行航向角 theta = %.4f rad = %.2f°\n', theta, rad2deg(theta));
fprintf('飞行速度   v     = %.4f m/s\n', v);
fprintf('投放时刻   t_l   = %.4f s\n', t_l);
fprintf('引信延时   tau   = %.4f s\n', tau);
fprintf('起爆时刻   t_det = %.4f s\n', t_det);
fprintf('投放点 P_l       = (%.2f, %.2f, %.2f) m\n', P_l);
fprintf('起爆点 P_d       = (%.2f, %.2f, %.2f) m\n', P_d);
fprintf('有效遮蔽时长     = %.4f s\n', T_shield);
fprintf('=============================================\n');

%% ================== 5. 遮蔽时间线图 ==================
t_full = 0 : params.dt : params.T_end;
shielded = false(size(t_full));
for k = 1:numel(t_full)
    t = t_full(k);
    if t >= t_det && t <= t_det + 20
        C = P_d - [0; 0; 3*(t - t_det)];
        M = params.M0 + params.v_m*t*params.u_m;
        shielded(k) = q2_occlusion(M, C, params.R, params);
    end
end

% 提取遮蔽区间
edges = diff([0, shielded, 0]);
starts = find(edges == 1) * params.dt - params.dt;
stops  = find(edges == -1) * params.dt - params.dt;
fprintf('遮蔽区间（共 %d 段）:\n', numel(starts));
for i = 1:numel(starts)
    fprintf('  [%.2f, %.2f] s，时长 %.2f s\n', starts(i), stops(i), stops(i)-starts(i));
end

figure('Color','w','Position',[80 80 820 420]);
hold on; grid on;
% 用 patch 填充遮蔽时段
for i = 1:numel(starts)
    fill([starts(i), stops(i), stops(i), starts(i)], ...
         [0, 0, 1, 1], [0.9 0.2 0.2], 'EdgeColor','none', ...
         'FaceAlpha', 0.25);
end
stairs(t_full, shielded, 'k-', 'LineWidth', 1.5);
yline(0,'-','未遮蔽');  yline(1,'-','遮蔽');
xline(t_det,      '--b', '起爆时刻');
xline(t_det + 20, '--b', '起爆+20s');
xline(params.T_end, '--r', '导弹到假目标');
xlabel('时间 t (s)');  ylabel('遮蔽状态');
ylim([-0.15 1.15]);  yticks([0 1]);
title(sprintf('M1 对真目标的遮蔽时间线（总遮蔽时长 = %.2f s）', T_shield));
hold off;

%% ================== 6. 三维态势图 ==================
figure('Color','w','Position',[100 100 900 680]);
hold on; grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('问题2 最优投放策略三维态势图');

% 导弹轨迹（红色）
mt = linspace(0, params.T_end, 300);
Mtraj = params.M0' + params.v_m*mt'*params.u_m';
plot3(Mtraj(:,1), Mtraj(:,2), Mtraj(:,3), 'r-', 'LineWidth', 1.6);
scatter3(params.M0(1), params.M0(2), params.M0(3), 90, 'r', '^', 'filled');
text(params.M0(1)+300, params.M0(2), params.M0(3), 'M1', 'Color', 'r', 'FontSize', 11);

% 无人机轨迹（蓝色）
ut = linspace(0, t_l, 100);
UAVtraj = params.uav0' + v*ut'*u_h';
plot3(UAVtraj(:,1), UAVtraj(:,2), UAVtraj(:,3), 'b-', 'LineWidth', 1.6);
scatter3(params.uav0(1), params.uav0(2), params.uav0(3), 80, 'b', 'o', 'filled');
text(params.uav0(1)+300, params.uav0(2), params.uav0(3), 'FY1', 'Color', 'b', 'FontSize', 11);

% 投放点与起爆点
scatter3(P_l(1), P_l(2), P_l(3), 100, 'g', 'd', 'filled');
scatter3(P_d(1), P_d(2), P_d(3), 100, 'k', 'p', 'filled');
text(P_l(1)+300, P_l(2), P_l(3), '投放点', 'Color', 'g', 'FontSize', 10);
text(P_d(1)+300, P_d(2), P_d(3), '起爆点', 'Color', 'k', 'FontSize', 10);

% 云团下沉轨迹（品红）
cdt = linspace(t_det, min(t_det+20, params.T_end), 80);
Ctraj = P_d - [0; 0; 3]*(cdt - t_det);
plot3(Ctraj(1,:), Ctraj(2,:), Ctraj(3,:), 'm-', 'LineWidth', 1.6);

% 假目标（原点）
scatter3(0, 0, 0, 90, 'm', '*');
text(500, 0, -100, '假目标(0,0,0)', 'Color', 'm', 'FontSize', 10);

% 真目标圆柱体（绿色半透明）
theta_c = linspace(0, 2*pi, 60);
x_cyl = params.cylCenter(1) + params.cylR*cos(theta_c);
y_cyl = params.cylCenter(2) + params.cylR*sin(theta_c);
Xc = repmat(x_cyl, 2, 1);
Yc = repmat(y_cyl, 2, 1);
Zc = [params.cylBottom*ones(size(x_cyl)); params.cylTop*ones(size(x_cyl))];
surf(Xc, Yc, Zc, 'FaceColor', 'g', 'FaceAlpha', 0.25, 'EdgeColor', 'g', 'LineWidth', 0.8);
text(params.cylCenter(1), params.cylCenter(2), params.cylTop+15, ...
     '真目标圆柱 r=7m h=10m', 'Color', 'g', 'FontSize', 10);

legend('导弹轨迹','M1','无人机轨迹','FY1','投放点','起爆点','云团下沉轨迹', ...
       '假目标','真目标','Location','best');
view(3);
hold off;

disp('求解完成。请将控制台输出的最优策略写入论文，并检查两张图。');
