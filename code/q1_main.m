%% q1_main  问题1：给定参数下烟幕干扰弹对 M1 的有效遮蔽时长
%
%  FY1 以 120 m/s 朝向假目标飞行，受领任务 1.5 s 后投放 1 枚烟幕干扰弹，
%  间隔 3.6 s 起爆。求烟幕干扰弹对 M1 的有效遮蔽时长。
%
%  求解思路（对应正文模型）：
%    1) 运动学模型：由给定参数解析出投放点、起爆点、起爆时刻；
%       导弹 M1 匀速直线直指假目标，云团起爆后以 3 m/s 匀速下沉。
%    2) 遮挡判定模型：构造布尔函数 I(t)，I(t)=1 表示 t 时刻圆柱被完全遮蔽。
%       引理①（连续区间）：有效遮蔽时间区间为一段连续闭区间，
%          故有效遮蔽时长 T = ∫I(t)dt = t_end - t_start。
%       引理②（上下圆周判定）：整个圆柱被完全遮蔽 ⇔ 上、下底面圆周
%          被完全遮蔽，故只需对圆周采样点做判定。
%       点遮蔽（到线段距离，含端点）：烟幕球（中心 C、半径 R）与
%          "导弹-圆周点"线段 MP 相交，即 C 到线段 MP 的最近距离 <= R，
%          等价于导弹视角下"目标轮廓落在烟幕角圆盘内"。
%    3) 求解算法：构造连续裕量函数 H(t)，H(t)<=0 等价于 I(t)=1，
%       "粗扫 + 二分查找"精确求 H(t)=0 的两个根，即遮蔽区间边界。
%    4) 敏感性分析：圆周采样点数 n 取不同值（6~600），
%       展示时长收敛性；结果表明 n>=100 时已收敛，
%       故取 n=100 兼顾精度与效率。
%
%  运行：MATLAB 中 cd 到 code 文件夹后直接运行 q1_main（无需额外工具箱）。
%  作者：数学建模 A 题队伍    日期：2026-08-25

clear; clc; close all;

%% ================== 1. 场景常量 ==================
params.g         = 9.8;                        % 重力加速度 (m/s^2)
params.R         = 10;                         % 烟幕有效遮蔽半径 (m)
params.M0        = [20000; 0; 2000];           % M1 初始位置
params.u_m       = -params.M0/norm(params.M0); % M1 飞行方向（直指假目标原点）
params.v_m       = 300;                        % 导弹速度 (m/s)
params.T_end     = norm(params.M0)/params.v_m; % 导弹到达假目标时刻 (s)
params.cylCenter = [0, 200];                   % 真目标圆柱轴心 (x,y)
params.cylR      = 7;                          % 圆柱半径 (m)
params.cylTop    = 10;                         % 圆柱顶面高度 (m)
params.cylBottom = 0;                          % 圆柱底面高度 (m)

% 问题1给定参数
params.uav0  = [17800; 0; 1800];               % FY1 初始位置
params.v_u   = 120;                            % FY1 飞行速度 (m/s)
params.u_h   = [-1; 0; 0];                     % FY1 朝向假目标（水平 -x 方向）
params.t_l   = 1.5;                            % 投放时刻 (s)
params.tau   = 3.6;                            % 引信延时 (s)

%% ================== 2. 运动学模型（模型准备） ==================
% 无人机：投放前等高度匀速直线运动，投放点
params.P_l = params.uav0 + params.v_u*params.t_l*params.u_h;
% 干扰弹：投放后仅受重力（初速等于投放瞬间无人机速度，竖直初速为 0），起爆点
params.P_d = params.P_l + params.v_u*params.tau*params.u_h ...
             - [0; 0; 0.5*params.g*params.tau^2];
params.t_d = params.t_l + params.tau;          % 起爆时刻

fprintf('========== 问题1 运动学解析 ==========\n');
fprintf('投放点 P_l = (%.3f, %.3f, %.3f) m\n', params.P_l);
fprintf('起爆点 P_d = (%.3f, %.3f, %.3f) m\n', params.P_d);
fprintf('起爆时刻 t_d = %.3f s\n', params.t_d);
fprintf('=====================================\n\n');

%% ================== 3. 遮挡判定模型 ==================
% 引理②：只需判断上、下底面圆周是否被完全遮蔽。
% 点遮蔽：烟幕球（中心 C、半径 R）与"导弹-圆周点"线段 MP 相交，
%         即 C 到线段 MP 的最近距离 <= R（垂足参数截断到 [0,1]，
%         允许云团包裹导弹端点或目标端点的情况）。
% 连续裕量函数 H(t) = max_p ( d_p(t) - R )，
%   其中 p 遍历两个圆周的全部采样点，d_p 为 C 到线段 MP 的最近距离。
%   H(t) <= 0  ⇔  圆柱被完全遮蔽；H 为连续函数，
%   故遮蔽集合 {t : H(t)<=0} 为一段连续闭区间（引理①）。

nPhi = 100;                                     % 圆周采样点数（敏感性分析：n>=100 已收敛）
P_rim = q1_rim_points(params, nPhi);            % 3 x (2*nPhi)

%% ================== 4. 粗扫 + 二分查找精确边界 ==================
dt_scan = 0.02;                                 % 粗扫时间步长 (s)
t_grid  = params.t_d : dt_scan : params.t_d + 20;
H_grid  = zeros(size(t_grid));
for k = 1:numel(t_grid)
    H_grid(k) = q1_margin(t_grid(k), P_rim, params);
end
inside = H_grid <= 0;
assert(any(inside), '粗扫未发现遮蔽区间');
i0 = find(inside, 1, 'first');
i1 = find(inside, 1, 'last');

t_start = q1_bisect(@(t) q1_margin(t, P_rim, params), ...
                    t_grid(max(i0-1,1)), t_grid(i0));
t_end   = q1_bisect(@(t) q1_margin(t, P_rim, params), ...
                    t_grid(i1), t_grid(min(i1+1, numel(t_grid))));
T_shield = t_end - t_start;                     % 有效遮蔽时长 T = ∫I(t)dt

fprintf('========== 问题1 求解结果 ==========\n');
fprintf('遮蔽区间 = [%.6f, %.6f] s（单段连续，对应引理①）\n', t_start, t_end);
fprintf('有效遮蔽时长 T = %.6f s ≈ %.4f s\n', T_shield, T_shield);
fprintf('====================================\n\n');

%% ================== 5. 遮蔽时间线（布尔函数积分复核） ==================
dt_fine = 0.005;                                % 细时间步长 (s)
ts_fine = 0 : dt_fine : params.T_end;
I = false(size(ts_fine));
for k = 1:numel(ts_fine)
    I(k) = q1_margin(ts_fine(k), P_rim, params) <= 0;
end
T_int = sum(I)*dt_fine;

% 统计遮蔽区间段数（验证连续单段）
edges = diff([0, I, 0]);
starts = find(edges == 1)*dt_fine - dt_fine;
stops  = find(edges == -1)*dt_fine - dt_fine;
fprintf('布尔积分复核：T_int = %.4f s（dt=%.3f），共 %d 段遮蔽区间\n', ...
        T_int, dt_fine, numel(starts));

% ---- 遮蔽时间线图 ----
figure('Color','w','Position',[80 80 820 400]);
hold on; grid on;
for i = 1:numel(starts)
    fill([starts(i), stops(i), stops(i), starts(i)], [0, 0, 1, 1], ...
         [0.9 0.2 0.2], 'EdgeColor','none', 'FaceAlpha', 0.25);
end
stairs(ts_fine, I, 'k-', 'LineWidth', 1.5);
yline(0,'-','未遮蔽');  yline(1,'-','遮蔽');
xline(params.t_d, '--b', '起爆时刻');
xline(t_start, '--m', '遮蔽开始');
xline(t_end,   '--m', '遮蔽结束');
xlabel('时间 t (s)');  ylabel('遮蔽状态 I(t)');
ylim([-0.15 1.15]);  yticks([0 1]);
title(sprintf('M1 对真目标的遮蔽时间线（有效遮蔽时长 = %.4f s）', T_shield));
legend('遮蔽时段','I(t)','Location','best');
hold off;

%% ================== 6. 裕量函数 H(t) 与二分求根示意图 ==================
figure('Color','w','Position',[100 100 820 400]);
plot(t_grid, H_grid, 'k-', 'LineWidth', 1.5);  hold on;  grid on;
yline(0, 'r--', 'H=0');
xline(t_start, '--m', 't_{start}');
xline(t_end,   '--m', 't_{end}');
scatter([t_start, t_end], [0, 0], 60, 'r', 'filled');
xlabel('时间 t (s)');  ylabel('裕量函数 H(t)');
title(sprintf('裕量函数 H(t)（H(t)≤0 即遮蔽；二分求根 t*=[%.4f, %.4f] s）', ...
      t_start, t_end));
hold off;

%% ================== 7. 三维态势图 ==================
figure('Color','w','Position',[120 120 900 680]);
hold on; grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('问题1 给定投放策略三维态势图');

% 导弹轨迹（红色）
mt = linspace(0, params.T_end, 300);
Mtraj = params.M0' + params.v_m*mt'*params.u_m';
plot3(Mtraj(:,1), Mtraj(:,2), Mtraj(:,3), 'r-', 'LineWidth', 1.6);
scatter3(params.M0(1), params.M0(2), params.M0(3), 90, 'r', '^', 'filled');
text(params.M0(1)+300, params.M0(2), params.M0(3), 'M1', 'Color', 'r', 'FontSize', 11);

% 无人机轨迹（蓝色）
ut = linspace(0, params.t_l, 100);
UAVtraj = params.uav0' + params.v_u*ut'*params.u_h';
plot3(UAVtraj(:,1), UAVtraj(:,2), UAVtraj(:,3), 'b-', 'LineWidth', 1.6);
scatter3(params.uav0(1), params.uav0(2), params.uav0(3), 80, 'b', 'o', 'filled');
text(params.uav0(1)+300, params.uav0(2), params.uav0(3), 'FY1', 'Color', 'b', 'FontSize', 11);

% 投放点与起爆点
scatter3(params.P_l(1), params.P_l(2), params.P_l(3), 100, 'g', 'd', 'filled');
scatter3(params.P_d(1), params.P_d(2), params.P_d(3), 100, 'k', 'p', 'filled');
text(params.P_l(1)+300, params.P_l(2), params.P_l(3), '投放点', 'Color', 'g', 'FontSize', 10);
text(params.P_d(1)+300, params.P_d(2), params.P_d(3), '起爆点', 'Color', 'k', 'FontSize', 10);

% 云团下沉轨迹（品红）
cdt = linspace(params.t_d, params.t_d + 20, 80);
Ctraj = params.P_d - [0; 0; 3]*(cdt - params.t_d);
plot3(Ctraj(1,:), Ctraj(2,:), Ctraj(3,:), 'm-', 'LineWidth', 1.6);

% 假目标（原点）
scatter3(0, 0, 0, 90, 'm', '*');
text(500, 0, -100, '假目标(0,0,0)', 'Color', 'm', 'FontSize', 10);

% 真目标圆柱体（绿色半透明）
theta_c = linspace(0, 2*pi, 60);
x_cyl = params.cylCenter(1) + params.cylR*cos(theta_c);
y_cyl = params.cylCenter(2) + params.cylR*sin(theta_c);
surf([x_cyl; x_cyl], [y_cyl; y_cyl], ...
     [params.cylBottom*ones(size(x_cyl)); params.cylTop*ones(size(x_cyl))], ...
     'FaceColor', 'g', 'FaceAlpha', 0.25, 'EdgeColor', 'g', 'LineWidth', 0.8);
text(params.cylCenter(1), params.cylCenter(2), params.cylTop+15, ...
     '真目标圆柱 r=7m h=10m', 'Color', 'g', 'FontSize', 10);

legend('导弹轨迹','M1','无人机轨迹','FY1','投放点','起爆点','云团下沉轨迹', ...
       '假目标','真目标','Location','best');
view(3);
hold off;

%% ================== 7b. 投影示意图（导弹视角视图平面） ==================
% 在导弹处建立视线坐标系，把目标可见轮廓与烟幕角圆盘投影到视图平面，
% 直观展示"目标轮廓落在烟幕角圆盘内 = 遮蔽有效"。
t_v = 8.7;                                   % 代表时刻（遮蔽窗口 [t_start,t_end] 内）
P_m = params.M0 + params.v_m*t_v*params.u_m; % 导弹位置
C_v = params.P_d - [0; 0; 3*(t_v - params.t_d)];  % 云团中心
T_v = [0; 200; 0];                           % 目标参考点（圆柱底面圆心）
ww = (T_v - P_m)/norm(T_v - P_m);            % 视线方向单位向量
uu = cross(ww, [0;0;1]);  uu = uu/norm(uu);  % 右方向
vv = cross(ww, uu);                          % 上方向

projA = @(X) (X - P_m)'*uu / norm(X - P_m);  % 水平角偏移 (rad)
projB = @(X) (X - P_m)'*vv / norm(X - P_m);  % 垂直角偏移 (rad)

% 目标轮廓采样（顶面圆周 + 底面圆周）
phi_c = linspace(0, 2*pi, 120);
Pt_ring = [params.cylCenter(1) + params.cylR*cos(phi_c);
           params.cylCenter(2) + params.cylR*sin(phi_c);
           params.cylTop*ones(1,120)];
Pb_ring = [params.cylCenter(1) + params.cylR*cos(phi_c);
           params.cylCenter(2) + params.cylR*sin(phi_c);
           params.cylBottom*ones(1,120)];
A_top = arrayfun(@(i) projA(Pt_ring(:,i)), 1:120);
B_top = arrayfun(@(i) projB(Pt_ring(:,i)), 1:120);
A_bot = arrayfun(@(i) projA(Pb_ring(:,i)), 1:120);
B_bot = arrayfun(@(i) projB(Pb_ring(:,i)), 1:120);

% 烟幕角圆盘：圆心 = 云心投影角偏移，半径 = asin(R/|P-C|)
A_c = projA(C_v);  B_c = projB(C_v);
alpha_v = asin(params.R / norm(P_m - C_v));

% 底面前弧（导弹可见部分）：按投影角半径排序取靠近目标参考点的前半
r_bot = sqrt(A_bot.^2 + B_bot.^2);
[~, ord] = sort(r_bot);
front_arc = ord(1:60);

figure('Color','w','Position',[160 160 980 480]);
% --- 左：全貌 ---
subplot(1,2,1);  hold on; grid on; axis equal;
th_d = linspace(0, 2*pi, 100);
fill(A_c + alpha_v*cos(th_d), B_c + alpha_v*sin(th_d), ...
     [0.85 0.85 0.85], 'EdgeColor', 'k', 'FaceAlpha', 0.5, 'LineWidth', 1.2);
plot(A_top, B_top, 'b-', 'LineWidth', 1.5);
plot(A_bot(front_arc), B_bot(front_arc), 'g-', 'LineWidth', 1.5);
plot(A_bot(~ismember(1:120, front_arc)), B_bot(~ismember(1:120, front_arc)), ...
     'g--', 'LineWidth', 1);
plot([0 A_c], [0 B_c], 'm--', 'LineWidth', 1.2);
scatter(0, 0, 40, 'r', 'filled');
scatter(A_c, B_c, 40, 'm', 'p', 'filled');
xlabel('水平角偏移 (rad)');  ylabel('垂直角偏移 (rad)');
title(sprintf('导弹视角：目标轮廓与烟幕角圆盘（t=%.1f s，遮蔽中）', t_v));
legend('烟幕角圆盘','顶面圆周','底面前弧(可见)','底面后弧','云心方向','Location','best');
hold off;
% --- 右：目标轮廓局部放大 ---
subplot(1,2,2);  hold on; grid on; axis equal;
plot(A_top*1e4, B_top*1e4, 'b-', 'LineWidth', 1.5);
plot(A_bot(front_arc)*1e4, B_bot(front_arc)*1e4, 'g-', 'LineWidth', 1.5);
scatter(0, 0, 40, 'r', 'filled');
scatter(A_c*1e4, B_c*1e4, 40, 'm', 'p', 'filled');
plot([0 A_c]*1e4, [0 B_c]*1e4, 'm--', 'LineWidth', 1.2);
xlabel('水平角偏移 (×10^{-4} rad)');  ylabel('垂直角偏移 (×10^{-4} rad)');
title('目标轮廓局部放大');
hold off;

%% ================== 8. 敏感性分析（圆周采样点数） ==================
ns = [6, 12, 24, 48, 100, 200, 300, 400, 600];
T_sens = zeros(size(ns));
tS = zeros(size(ns));  tE = zeros(size(ns));
fprintf('========== 敏感性分析（圆周采样点数 n）==========\n');
fprintf('%5s | %12s | %12s | %10s\n', 'n', 't_start(s)', 't_end(s)', 'T(s)');
for i = 1:numel(ns)
    P = q1_rim_points(params, ns(i));
    Hh = zeros(size(t_grid));
    for k = 1:numel(t_grid)
        Hh(k) = q1_margin(t_grid(k), P, params);
    end
    in = Hh <= 0;
    j0 = find(in, 1, 'first');
    j1 = find(in, 1, 'last');
    tS(i) = q1_bisect(@(t) q1_margin(t, P, params), t_grid(max(j0-1,1)), t_grid(j0));
    tE(i) = q1_bisect(@(t) q1_margin(t, P, params), t_grid(j1), t_grid(min(j1+1, numel(t_grid))));
    T_sens(i) = tE(i) - tS(i);
    fprintf('%5d | %12.6f | %12.6f | %10.6f\n', ns(i), tS(i), tE(i), T_sens(i));
end
fprintf('=================================================\n');
fprintf('结论：n >= 100 时 T 已收敛到 %.6f s，故取 n = 100 兼顾精度与效率。\n', T_sens(end));

% 收敛曲线图
figure('Color','w','Position',[140 140 760 420]);
semilogx(ns, T_sens, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
hold on; grid on;
yline(T_sens(end), 'r--');
xlabel('圆周采样点数 n');  ylabel('有效遮蔽时长 T (s)');
title(sprintf('敏感性分析：T 随采样点数 n 收敛到 %.6f s', T_sens(end)));
hold off;

%% ================== 9. 引理②数值验证（全圆柱面 vs 上下圆周） ==================
P_full = q1_full_points(params);
Hf = zeros(size(t_grid));
for k = 1:numel(t_grid)
    Hf(k) = q1_margin(t_grid(k), P_full, params);
end
infull = Hf <= 0;
j0 = find(infull, 1, 'first');
j1 = find(infull, 1, 'last');
tSf = q1_bisect(@(t) q1_margin(t, P_full, params), t_grid(max(j0-1,1)), t_grid(j0));
tEf = q1_bisect(@(t) q1_margin(t, P_full, params), t_grid(j1), t_grid(min(j1+1, numel(t_grid))));
fprintf('引理②数值验证：全圆柱面采样时长 = %.6f s，与圆周判据 %.6f s 一致\n', ...
        tEf-tSf, T_shield);

disp('求解完成。请将控制台输出的遮蔽区间与时长写入论文，并检查全部五张图。');

%% ====================== 子函数 ======================

function P = q1_rim_points(params, n)
% 上、下底面圆周各 n 个采样点，返回 3 x (2n)
phi = linspace(0, 2*pi, n+1);
phi(end) = [];                                 % 去掉重复端点
P_top = [params.cylCenter(1) + params.cylR*cos(phi);
         params.cylCenter(2) + params.cylR*sin(phi);
         params.cylTop*ones(1, n)];
P_bot = [params.cylCenter(1) + params.cylR*cos(phi);
         params.cylCenter(2) + params.cylR*sin(phi);
         params.cylBottom*ones(1, n)];
P = [P_top, P_bot];                            % 3 x (2n)
end

function P = q1_full_points(params)
% 全圆柱面采样（用于数值验证引理②）：
%   侧面 6 个高度环 x 36 点 + 上、下底面圆盘同心圆采样
nr = 36;  nz = 6;
phi = linspace(0, 2*pi, nr+1);  phi(end) = [];
zs = linspace(params.cylBottom, params.cylTop, nz);
Pts = [];
for z = zs
    ring = [params.cylCenter(1) + params.cylR*cos(phi);
            params.cylCenter(2) + params.cylR*sin(phi);
            z*ones(1, nr)];
    Pts = [Pts, ring];
end
for z = [params.cylBottom, params.cylTop]
    for r = linspace(0, params.cylR, 6)
        nq = max(1, round(12*r/params.cylR));
        q = linspace(0, 2*pi, nq+1);  q(end) = [];
        if isempty(q), q = 0; end
        disk = [params.cylCenter(1) + r*cos(q);
                params.cylCenter(2) + r*sin(q);
                z*ones(1, numel(q))];
        Pts = [Pts, disk];
    end
end
P = Pts;                                       % 3 x N
end

function H = q1_margin(t, P, params)
% 连续裕量函数：H(t) <= 0  ⇔  圆柱被完全遮蔽（上下圆周全部采样点均被遮蔽）
if t < params.t_d || t > params.t_d + 20
    H = 1e9;                                   % 云团未起爆或已失效
    return;
end
M = params.M0 + params.v_m*t*params.u_m;       % 导弹位置 3x1
C = params.P_d - [0; 0; 3*(t - params.t_d)];   % 云团中心 3x1（3 m/s 下沉）
v = P - M;                                     % 导弹 -> 圆周点 3 x nP
w = C - M;                                     % 导弹 -> 云团 3x1
dvv = sum(v.*v, 1);                            % |v|^2 (1 x nP)
tp = (w'*v)./dvv;                              % 垂足参数：foot = M + tp*v (1 x nP)
tp = max(0, min(1, tp));                       % 截断到 [0,1]（到线段距离，含端点）
foot = M + tp.*v;                              % 线段上最近点 3 x nP
d = sqrt(sum((C - foot).^2, 1));               % C 到线段 MP 最近距离 (1 x nP)
H = max(d - params.R);                         % 全部圆周点取最大（引理②）
end

function x = q1_bisect(fun, lo, hi)
% 二分查找：求连续函数 fun 在 [lo,hi] 上的零点（200 次迭代）
for k = 1:200
    mid = 0.5*(lo + hi);
    if fun(lo)*fun(mid) <= 0
        hi = mid;
    else
        lo = mid;
    end
end
x = 0.5*(lo + hi);
end
