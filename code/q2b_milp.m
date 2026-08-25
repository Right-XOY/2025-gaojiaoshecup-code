function q2b_milp
% q2b_milp  问题2 方案B：MILP 离散网格 + 局部精细搜索 + 精确区间 + 全圆柱保守校验
%
% 对应 Python 版 code/q2b_milp.py，改用 MATLAB 内置 intlinprog
% （Optimization Toolbox，无需 Gurobi / HiGHS）。
%
% 决策变量（离散网格）：
%   FY1 航向角 theta、飞行速度 v、投放时刻 t_drop（投放后 3.6 s 起爆），
%   起爆后形成半径 10 m 球体，以 3 m/s 匀速下沉，起爆后 20 s 内有效。
%
% MILP 结构：
%   z_k : 0-1，是否选择第 k 个候选投放方案（恰好选 1 个）
%   s_j : 0-1，离散时刻 j 是否有效遮蔽
%   s_j <= sum_{k 能在 j 时刻遮蔽} z_k
%   目标 : max sum(s_j)，即最大化有效遮蔽时间
%
% 遮蔽判定（MILP 内）：云团中心到 "M1(t)-真目标底面圆心(0,200,0)" 视线段的
% 距离 <= 10 m（到线段距离含端点，与主判据口径一致）。为控制 MILP 规模，
% 用目标参考点（底面圆心）单点判定加速，求解后另做两件事：
%   1) 在最优网格解附近做局部精细搜索（航向/速度/投放时刻）；
%   2) 用上下圆周各 100 点校验 "整个圆柱完全遮蔽" 的时长（n>=100 收敛口径）。
%
% 运行：MATLAB 中 cd 到 code 文件夹后运行 q2b_milp（需 Optimization Toolbox）
% 作者：数学建模 A 题队伍    日期：2026-08-25

clear; clc;
rng(2026);

%% ===================== 1. 问题参数 =====================
G = 9.8;
M0 = [20000, 0, 2000];            % M1 初始位置
TGT = [0, 200, 0];                % 真目标底面圆心
F0 = [17800, 0, 1800];            % FY1 初始位置
V_M = -300 * M0 / norm(M0);       % M1 速度（指向假目标原点）
V_MIN = 70;  V_MAX = 140;         % FY1 速度范围
T_BOMB = 3.6;                     % 投放到起爆的间隔
R_SMOKE = 10.0;                   % 烟幕有效半径
V_SINK = 3.0;                     % 云团下沉速度
T_LIFE = 20.0;                    % 云团有效时长
DZ_DET = 0.5 * G * T_BOMB^2;      % 起爆前自由落体下降量
T_END = 67.0;                     % M1 到达原点时刻

HEAD_STEP = 6.0;                  % 航向角步长（度）
V_STEP = 10.0;                    % 速度步长
T_DROP_STEP = 0.5;                % 投放时刻步长
T_DROP_MAX = 20.0;
DT = 0.25;                        % MILP 时间离散步长

headings = 0 : HEAD_STEP : 360 - 1e-9;
speeds   = V_MIN : V_STEP : V_MAX;
t_drop_grid = 0 : T_DROP_STEP : T_DROP_MAX;
time_grid   = 0 : DT : T_END - 1e-9;

N_H = numel(headings);
N_V = numel(speeds);
N_R = numel(t_drop_grid);
K = N_H * N_V * N_R;
N_TIME = numel(time_grid);
fprintf('候选方案数 K = %d x %d x %d = %d，离散时刻数 = %d\n', ...
        N_H, N_V, N_R, K, N_TIME);

%% ===================== 2. 预计算 blockers =====================
% 展平顺序与 Python 一致：k = v*(N_H*N_R) + h*N_R + r（r 最快，列优先展平）
[Rg, Hg, Vg] = ndgrid(0:N_R-1, 0:N_H-1, 0:N_V-1);
r_lin = reshape(Rg, 1, []);                     % 1 x K 行向量（r 最快）
h_lin = reshape(Hg, 1, []);
v_lin = reshape(Vg, 1, []);
theta_k = headings(h_lin + 1);                  % 1 x K 航向角
v_k     = speeds(v_lin + 1);                    % 1 x K 速度
td_k    = t_drop_grid(r_lin + 1);               % 1 x K 投放时刻
vu = [v_k .* cosd(theta_k);                     % 3 x K 无人机水平速度矢量
      v_k .* sind(theta_k);
      zeros(1, K)];
det_t   = td_k + T_BOMB;                % 起爆时刻 1 x K
det_pos = F0(:) * ones(1, K) + vu .* det_t + [0; 0; -DZ_DET];   % 起爆点 3 x K

blockers = cell(N_TIME, 1);
for j = 1:N_TIME
    tj = time_grid(j);
    mj = M0(:) + V_M(:) * tj;                        % 导弹位置 3x1
    sink = zeros(3, K);                              % 云团下沉量 3xK
    sink(3, :) = -V_SINK * max(tj - det_t, 0);
    cj = det_pos + sink;                             % 云团中心 3xK
    active = (tj >= det_t) & (tj <= det_t + T_LIFE); % 有效期内 1xK
    d = dist_to_segment(cj, mj, TGT(:));             % 云团到视线段距离 1xK
    block = (d <= R_SMOKE) & active;
    blockers{j} = find(block);
end
fprintf('blockers 预计算完成\n');

%% ===================== 3. MILP 求解（intlinprog） =====================
% 决策变量 x = [z(1..K); s(1..N_TIME)]，全部 0-1
% 目标 max sum(s_j) => min f'x, f = [0; -1]
f = [zeros(K, 1); -ones(N_TIME, 1)];
intcon = 1 : (K + N_TIME);
lb = zeros(K + N_TIME, 1);
ub = ones(K + N_TIME, 1);

% 无候选遮蔽的时刻强制 s_j = 0（上界设为 0）
for j = 1:N_TIME
    if isempty(blockers{j})
        ub(K + j) = 0;
    end
end

% 等式约束：恰好选择 1 个候选方案  sum(z) = 1
Aeq = [ones(1, K), zeros(1, N_TIME)];
beq = 1;

% 不等式约束：s_j - sum_{k in blockers{j}} z_k <= 0
rows = [];  cols = [];  vals = [];
for j = 1:N_TIME
    idx = blockers{j};
    if isempty(idx), continue; end
    rows = [rows, j * ones(1, numel(idx) + 1)];
    cols = [cols, idx, K + j];
    vals = [vals, -ones(1, numel(idx)), 1];
end
A = sparse(rows, cols, vals, N_TIME, K + N_TIME);
b = zeros(N_TIME, 1);

opts = optimoptions('intlinprog', 'Display', 'off', 'MaxTime', 120);
[xopt, ~, exitflag] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, opts);
if exitflag <= 0
    warning('intlinprog 未求出最优解，exitflag = %d', exitflag);
end

best_k = find(xopt(1:K) > 0.5, 1) - 1;       % 转回 0 基下标
obj_samples = round(sum(xopt(K+1:end)));     % 有效遮蔽离散时刻数
fprintf('求解器: intlinprog，MILP 目标值 = %d 个离散时刻，约 %.2f s\n', ...
        obj_samples, obj_samples * DT);

%% ===================== 4. 提取最优网格方案 =====================
r_idx = mod(best_k, N_R);
tmp   = floor(best_k / N_R);
h_idx = mod(tmp, N_H);
v_idx = floor(tmp / N_H);

theta_star = headings(h_idx + 1);
v_star     = speeds(v_idx + 1);
td_star    = t_drop_grid(r_idx + 1);
[release0, det_t0, det0] = candidate_geometry(theta_star, v_star, td_star, ...
    F0, T_BOMB, DZ_DET);

fprintf('\n===== MILP 网格最优方案 =====\n');
fprintf('航向角 theta = %.1f deg，速度 v = %.1f m/s\n', theta_star, v_star);
fprintf('投放时刻 t_drop = %.2f s\n', td_star);
fprintf('投放点 = (%.2f, %.2f, %.2f)\n', release0(1), release0(2), release0(3));
fprintf('起爆时刻 t_det = %.2f s\n', det_t0);
fprintf('起爆点 = (%.2f, %.2f, %.2f)\n', det0(1), det0(2), det0(3));

%% ===================== 5. 局部精细搜索 =====================
fprintf('\n===== 局部精细搜索 =====\n');
best = [theta_star, v_star, td_star];
best_dur = blocked_duration(best(1), best(2), best(3), ...
    M0, V_M, TGT, F0, T_BOMB, DZ_DET, V_SINK, R_SMOKE, T_LIFE, T_END);

for theta = theta_star - 3 : 0.5 : theta_star + 3
    for v = max(V_MIN, v_star - 10) : 1 : min(V_MAX, v_star + 10)
        for td = max(0, td_star - 0.5) : 0.05 : td_star + 0.5
            dur = blocked_duration(theta, v, td, ...
                M0, V_M, TGT, F0, T_BOMB, DZ_DET, V_SINK, R_SMOKE, T_LIFE, T_END);
            if dur > best_dur
                best_dur = dur;
                best = [theta, v, td];
            end
        end
    end
end
theta_s = best(1);  v_s = best(2);  td_s = best(3);
[release_s, det_t_s, det_s] = candidate_geometry(theta_s, v_s, td_s, ...
    F0, T_BOMB, DZ_DET);
fprintf('精细最优：theta = %.2f deg，v = %.2f m/s，t_drop = %.2f s\n', ...
        theta_s, v_s, td_s);
fprintf('精细时长（目标参考点单点判据，步长 0.02 s 粗测）= %.3f s\n', best_dur);

%% ===================== 6. 二分法精确求遮蔽区间（目标参考点单点判据） =====================
dcenter = @(t) dcenter_func(t, det_t_s, det_s, T_LIFE, V_SINK, ...
    M0, V_M, TGT, R_SMOKE);

tscan = max(0, det_t_s) : 0.01 : min(T_END, det_t_s + T_LIFE);
dscan = zeros(size(tscan));
for i = 1:numel(tscan)
    dscan(i) = dcenter(tscan(i));
end
inside = dscan <= R_SMOKE;
t_start_exact = NaN;  t_end_exact = NaN;
if any(inside)
    i0 = find(inside, 1, 'first');
    i1 = find(inside, 1, 'last');
    t_start_exact = bisect(tscan(max(1, i0-1)), tscan(i0), dcenter, R_SMOKE);
    t_end_exact   = bisect(tscan(i1), tscan(min(i1+1, numel(tscan))), dcenter, R_SMOKE);
end

if ~isnan(t_start_exact)
    dur_exact = t_end_exact - t_start_exact;
    fprintf('\n===== 精确有效遮蔽（目标参考点单点判据） =====\n');
    fprintf('遮蔽区间: [%.4f, %.4f] s\n', t_start_exact, t_end_exact);
    fprintf('有效遮蔽时长: %.4f s\n', dur_exact);
else
    dur_exact = 0;
    fprintf('目标参考点单点判据下无有效遮蔽\n');
end

%% ===================== 7. 全圆柱完全遮蔽校验（上下圆周各 100 点，n>=100 收敛口径） =====================
N_RIM = 100;
ang = (0:N_RIM-1) * 2 * pi / N_RIM;              % 100 个均匀点（不含端点，与 Python endpoint=False 一致）
rim = zeros(2 * N_RIM, 3);
k = 0;
for zz = [0, 10]
    for a = ang
        k = k + 1;
        rim(k, :) = [7 * cos(a), 200 + 7 * sin(a), zz];
    end
end
rim = rim';                                      % 3 x (2*N_RIM)

ts_rim = max(0, det_t_s) : 0.02 : min(T_END, det_t_s + T_LIFE);
mask_rim = false(1, numel(ts_rim));
for i = 1:numel(ts_rim)
    t = ts_rim(i);
    m = M0(:) + V_M(:) * t;
    c = det_s(:) + [0; 0; -V_SINK * (t - det_t_s)];
    d = dist_to_segment(rim, c, m);
    mask_rim(i) = max(d) <= R_SMOKE;
end
dur_rim = sum(mask_rim) * 0.02;
fprintf('\n===== 全圆柱完全遮蔽校验（上下圆周各 %d 点） =====\n', N_RIM);
fprintf('完全遮蔽时长 ≈ %.3f s\n', dur_rim);

%% ===================== 8. 汇总输出 =====================
fprintf('\n===== 最终策略汇总 =====\n');
fprintf('FY1 航向: 自 x 轴正向逆时针 %.2f°（即向假目标方向偏转 %.2f°）\n', ...
        theta_s, 180 - theta_s);
fprintf('FY1 速度: %.2f m/s\n', v_s);
fprintf('投放点: (%.2f, %.2f, %.2f)，t = %.2f s\n', ...
        release_s(1), release_s(2), release_s(3), td_s);
fprintf('起爆点: (%.2f, %.2f, %.2f)，t = %.2f s\n', ...
        det_s(1), det_s(2), det_s(3), det_t_s);
fprintf('有效遮蔽时长（新遮蔽判定：上下圆周完全遮蔽）约 %.3f s\n', dur_rim);
fprintf('（MILP 目标参考点单点判据时长 %.4f s，仅作求解加速参考）\n', dur_exact);

end

%% ===================== 子函数 =====================

function d = dist_to_segment(p, a, b)
% 向量化：p 为 3xK，a/b 为 3x1，返回到线段 ab 的距离 1xK
ab = b - a;                                                    % 3x1（或 3xN）
t = max(0, min(1, sum((p - a) .* ab, 1) ./ sum(ab .* ab, 1))); % 1xK 垂足参数裁剪
closest = a + t .* ab;                                         % 3xK（a 广播）
d = sqrt(sum((p - closest).^2, 1));                            % 1xK
end

function [release, det_t, det] = candidate_geometry(theta, v, t_drop, F0, T_BOMB, DZ_DET)
% 由 (theta, v, t_drop) 计算投放点、起爆时刻、起爆点
vu = v * [cosd(theta); sind(theta); 0];
release = F0(:) + vu * t_drop;
det_t = t_drop + T_BOMB;
det = release + vu * T_BOMB + [0; 0; -DZ_DET];
end

function dur = blocked_duration(theta, v, t_drop, M0, V_M, TGT, F0, T_BOMB, ...
    DZ_DET, V_SINK, R_SMOKE, T_LIFE, T_END)
% 用 0.02 s 细步长计算目标参考点单点判据下的有效遮蔽时长（MILP 求解加速代理）
[~, det_t, det] = candidate_geometry(theta, v, t_drop, F0, T_BOMB, DZ_DET);
t1 = max(0, det_t);
t2 = min(T_END, det_t + T_LIFE);
if t2 <= t1, dur = 0; return; end
ts = t1 : 0.02 : t2;
if isempty(ts), dur = 0; return; end
ms = M0(:) + V_M(:) * ts;                                % 3xN 导弹轨迹
sink = zeros(3, numel(ts));
sink(3, :) = -V_SINK * (ts - det_t);
cs = det(:) + sink;                                      % 3xN 云团中心
d = dist_to_segment(cs, ms, TGT(:));
dur = sum(d <= R_SMOKE) * 0.02;
end

function val = dcenter_func(t, det_t_s, det_s, T_LIFE, V_SINK, M0, V_M, TGT, R_SMOKE)
% 云团中心到 "M1(t)-真目标底面圆心" 视线段的距离；有效窗外返回大数
if t < det_t_s || t > det_t_s + T_LIFE
    val = 1e9;
    return;
end
m = M0(:) + V_M(:) * t;
c = det_s(:) + [0; 0; -V_SINK * (t - det_t_s)];
val = dist_to_segment(c, m, TGT(:));
end

function xr = bisect(lo, hi, f, R_SMOKE)
% 二分求 f(t) - R_SMOKE = 0 的根（f 返回距离）
for k = 1:200
    mid = 0.5 * (lo + hi);
    if (f(lo) - R_SMOKE) * (f(mid) - R_SMOKE) <= 0
        hi = mid;
    else
        lo = mid;
    end
end
xr = 0.5 * (lo + hi);
end
