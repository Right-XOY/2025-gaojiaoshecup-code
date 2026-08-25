%% q5_main  问题5：5 架无人机、每机至多 3 弹，干扰 M1/M2/M3 三枚来袭导弹
%  给出各无人机的飞行方向、飞行速度，各烟幕干扰弹的投放点与起爆点，
%  使三枚导弹各自的遮蔽时长之和（并集主指标）尽可能长，
%  并将结果写入 data/result3.xlsx。
%
%  求解架构：分层优化（四层）
%    第0层【精确预扫描】code/q5_prescan.py 对 5 机 x 3 导弹的每个组合
%      做网格扫描，生成遮蔽潜力矩阵 pot 与最优单弹参数表（不硬编码）。
%    第1层【任务指派】基于 pot 矩阵自动求解整数线性规划
%      （ILP 框架，c_ij 取精确预扫描 pot），把 5 机的
%      弹位指派给 M1/M2/M3，取最优 + 次优共 k 套候选，全部向下传递；
%    第2层【子问题优化】对每套候选，M1/M2/M3 各自独立求解一个
%      多弹接力子问题（最大化该导弹独立遮蔽并集时长，为第3层提供
%      遮蔽窗口种子）；
%    第3层【联合精化】固定每机航向/速度（第2层结果），对全部弹位的
%      投放/延时做联合优化，目标 = 三导弹遮蔽时长之和（并集主指标）；
%    第4层【细步长精化·Memetic】以最优候选为基础，细步长小种群 JADE
%      在保持 θ/v 不变下继续微调投放/延时，拉长遮蔽时长之和。
%  最后：对每套候选评估"三导弹遮蔽时长之和"，取最大者为最终方案。
%
%  依赖函数：q5_assign.m、q5_missile_fitness.m、q5_joint_fitness.m、
%            q3_jade.m、q3_occlusion_window.m
%  运行：MATLAB 中 cd 到 code 文件夹后运行 q5_main（无需额外工具箱）
%
%  作者：数学建模 A 题队伍    日期：2026-08-24

clear; clc; close all;
rng(2026);                      % 固定随机种子，保证可复现

%% ================== 1. 场景常量 ==================
params.g         = 9.8;
params.R         = 10;
params.v_min     = 70;
params.v_max     = 140;
params.v_m       = 300;
params.cylCenter = [0, 200];
params.cylR      = 7;
params.cylTop    = 10;
params.cylBottom = 0;
params.nPhi      = 100;                            % 圆周采样点数（敏感性分析：n>=100 收敛）

% 三枚导弹初始位置（3x3，每列一枚：M1、M2、M3）
params.M0_all = [20000  19000  18000;
                    0     600   -600;
                 2000    2100   1900];
params.u_m_all = zeros(3, 3);
for m = 1:3
    params.u_m_all(:, m) = -params.M0_all(:, m)/norm(params.M0_all(:, m));
end
params.T_end = zeros(1, 3);
for m = 1:3
    params.T_end(m) = norm(params.M0_all(:, m))/params.v_m;
end

% 五架无人机初始位置（3x5，每列一架：FY1~FY5）
params.uav0_all = [17800  12000   6000  11000  13000;
                      0    1400  -3000   2000  -2000;
                   1800    1400    700   1800   1300];
uav_names = {'FY1', 'FY2', 'FY3', 'FY4', 'FY5'};
m_names   = {'M1', 'M2', 'M3'};

%% ================== 2. 遮蔽潜力矩阵（预扫描数据文件，不硬编码） ==================
% pot(uav, m) = 第 uav 架无人机单弹对第 m 枚导弹的最优遮蔽时长 (s)，
% 由 code/q5_prescan.py 精确网格扫描生成（θ 每 10°、v 每 10 m/s、
% t_l/τ 细网格），写入 data/q5_prescan_pot.csv，此处读取。
data_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
pot = readmatrix(fullfile(data_dir, 'q5_prescan_pot.csv'));

% 精确预扫描最优单弹参数（uav, m, theta, v, t_l, tau, bestT, win_s, win_e），
% 由 q5_prescan.py 生成（data/q5_prescan_par.csv），作为接力预扫描/第2层种子。
scan_par = readmatrix(fullfile(data_dir, 'q5_prescan_par.csv'));

params.dt = 0.05;                      % 粗时间步长 (s)

%% ================== 3. 接力预扫描 + 第1层：候选任务指派 ==================
% 对每个可行组合求「单机 3 弹接力」最优时长，作为 ILP 系数，修正旧假设
% 「同机多弹无法接力」——单机接力可显著延长遮蔽（FY2→M2 3.9→7.5s），
% 避免多余机（如 FY4）被错误指派到已被单机接力覆盖的导弹上造成浪费。
[pot_relay, par_relay] = q5_relay_prescan(params, pot, scan_par);

k = 3;
cand = q5_assign(pot_relay, k);

%% ================== 4. 第2层+第3层：逐候选求解 ==================
opts_s.NP = 30;  opts_s.MAXITER = 200; opts_s.p = 0.05; opts_s.c = 0.1;
opts_j.NP = 40;  opts_j.MAXITER = 250; opts_j.p = 0.05; opts_j.c = 0.1;

% 结果容器
best_info = struct('total', -inf, 'x_lower', [], 'hist_lower', [], ...
                   'x_joint', [], 'f_joint', [], 'ginfo', [], 'theta_v_all', []);

for t = 1:k
    fprintf('\n================ 候选指派 %d/%d ================\n', t, k);

    % ---- 展开候选：各导弹服务表 -> servers 结构 + 全局弹位表 ----
    servers = cell(1, 3);               % servers{m}
    ginfo.slots = zeros(0, 3);          % [导弹, 无人机, 弹位编号]
    for m = 1:3
        srv = cand(t).missile{m};
        servers{m}.uavs  = unique(srv(:, 1))';
        slots = [];
        for r = 1:size(srv, 1)
            uav_idx = srv(r, 1);  nb = srv(r, 2);
            for q = 1:nb
                slots(end+1, :) = [uav_idx, q];
                ginfo.slots(end+1, :) = [m, uav_idx, q];
            end
        end
        servers{m}.slots = slots;
    end

    % ---- 第2层：三枚导弹独立子问题 JADE ----
    % 每机 θ/v 纳入决策变量（同机多弹共享），与各弹位 t_l/τ 联合搜索，
    % 使 JADE 在单弹最优方向（迎向导弹，即接力方向）上优化速度/延时，
    % 实现同机多弹真正接力（满足 θ/v 一旦确定不再调整的题目约束）。
    x_lower  = cell(1, 3);
    f_lower  = zeros(1, 3);
    hist_lower = cell(1, 3);
    theta_v_all = zeros(5, 2);          % 每机 [theta, v]（第2层优化结果）

    for m = 1:3
        nU = numel(servers{m}.uavs);
        nS = size(servers{m}.slots, 1);
        if nS == 0, x_lower{m} = []; f_lower(m) = 0; continue; end

        % 决策变量 = [θ/v (2*nU), t_l/τ (2*nS)]
        lb = [];  ub = [];
        for u = 1:nU
            lb = [lb, 0, params.v_min];        % θ ∈ [0,360]
            ub = [ub, 360, params.v_max];      % v ∈ [70,140]
        end
        for s = 1:nS
            uav_idx = servers{m}.slots(s, 1);
            tau_max = sqrt(2*params.uav0_all(3, uav_idx)/params.g);
            lb = [lb, 0, 0];
            ub = [ub, params.T_end(m), tau_max];
        end

        % 种子：预扫描单弹 θ/v + 接力速度档 + 接力预扫描最优解（par_relay）
        seeds = q5_build_seeds(m, servers{m}, scan_par, params, par_relay);
        opts_s.X0 = seeds;
        fun_m = @(x) q5_missile_fitness(x, params, m, servers{m});
        D = 2*nU + 2*nS;
        fprintf('  M%d 子问题（%d 机 %d 弹，%d 维 JADE）...\n', m, nU, nS, D);
        [x_lower{m}, f_lower(m), hist_lower{m}] = q3_jade(fun_m, lb, ub, opts_s);

        % 从第2层解提取每机 θ/v，写入全局 theta_v_all
        for u = 1:nU
            uav_idx = servers{m}.uavs(u);
            theta_v_all(uav_idx, :) = [x_lower{m}(2*u-1), x_lower{m}(2*u)];
        end
    end
    total_lower = sum(-f_lower);
    fprintf('  第2层结果：三导弹遮蔽时长之和（子问题参考）= %.4f s (M1=%.2f M2=%.2f M3=%.2f)\n', ...
            total_lower, -f_lower(1), -f_lower(2), -f_lower(3));

    % ---- 第3层：联合精化（固定 theta/v，只优化各弹位 t_l/tau）----
    nS_all = size(ginfo.slots, 1);
    lb_j = [];  ub_j = [];
    x0_j = zeros(1, 2*nS_all);
    for s = 1:nS_all
        m_idx = ginfo.slots(s, 1);
        uav_idx = ginfo.slots(s, 2);
        tau_max = sqrt(2*params.uav0_all(3, uav_idx)/params.g);
        lb_j = [lb_j, 0, 0];
        ub_j = [ub_j, params.T_end(m_idx), tau_max];
    end
    % 初始种子：从第2层各子问题解中提取对应弹位的 t_l/tau
    for m = 1:3
        if isempty(x_lower{m}), continue; end
        nU_m = numel(servers{m}.uavs);
        nS_m = size(servers{m}.slots, 1);
        for s = 1:nS_m
            uav_idx = servers{m}.slots(s, 1);
            t_l = x_lower{m}(2*nU_m + 2*s - 1);
            tau = x_lower{m}(2*nU_m + 2*s);
            row = find(ginfo.slots(:, 1) == m & ginfo.slots(:, 2) == uav_idx ...
                     & ginfo.slots(:, 3) == servers{m}.slots(s, 2));
            x0_j(2*row - 1) = t_l;   x0_j(2*row) = tau;
        end
    end
    % 种子矩阵：x0 及其邻域变体
    seeds_j = repmat(x0_j, 5, 1);
    for r = 2:5
        seeds_j(r, 1:2:end) = seeds_j(r, 1:2:end) + (r-2)*0.4;   % t_l 微移
        seeds_j(r, 2:2:end) = seeds_j(r, 2:2:end) + (r-2)*0.2;   % tau 微移
    end
    opts_j.X0 = seeds_j;
    fun_j = @(x) q5_joint_fitness(x, params, ginfo, theta_v_all);
    fprintf('  第3层联合精化（%d 维 JADE）...\n', 2*nS_all);
    [x_joint, f_joint, hist_joint] = q3_jade(fun_j, lb_j, ub_j, opts_j);
    total_joint = -f_joint;
    fprintf('  第3层结果：三导弹遮蔽时长之和 = %.4f s\n', total_joint);

    % ---- 记录当前候选最优信息 ----
    if total_joint > best_info.total
        best_info.total       = total_joint;
        best_info.x_lower     = x_lower;
        best_info.hist_lower  = hist_lower;
        best_info.x_joint     = x_joint;
        best_info.f_joint     = f_joint;
        best_info.ginfo       = ginfo;
        best_info.theta_v_all = theta_v_all;
        best_info.servers     = servers;
        best_info.total_lower = total_lower;
    end
end

fprintf('\n========== 最优候选：三导弹遮蔽时长之和 = %.4f s ==========\n', ...
        best_info.total);

%% ================== 5. 第4层：细步长精化（Memetic 局部微调） ==================
params.dt = 0.005;                     % 精细时间步长 (s)
x_final = best_info.x_joint;
ginfo   = best_info.ginfo;
theta_v_all = best_info.theta_v_all;

opts_f.NP = 15;  opts_f.MAXITER = 120; opts_f.p = 0.05; opts_f.c = 0.1;
opts_f.X0 = repmat(x_final, 4, 1);     % 以第3层结果为种子，保持 θ/v 不变
for r = 2:4
    opts_f.X0(r, 1:2:end) = opts_f.X0(r, 1:2:end) + (r-1)*0.05;
    opts_f.X0(r, 2:2:end) = opts_f.X0(r, 2:2:end) + (r-1)*0.03;
end
fun_fine = @(x) q5_joint_fitness(x, params, ginfo, theta_v_all);
fprintf('========== 第4层：细步长精化中 ... ==========\n');
% 边界按各弹位所属导弹/无人机的物理上限构造
lb_f = zeros(1, 2*size(ginfo.slots, 1));
ub_f = zeros(1, 2*size(ginfo.slots, 1));
for s = 1:size(ginfo.slots, 1)
    m_idx = ginfo.slots(s, 1);  uav_idx = ginfo.slots(s, 2);
    ub_f(2*s-1) = params.T_end(m_idx);
    ub_f(2*s)   = sqrt(2*params.uav0_all(3, uav_idx)/params.g);
end
[x_final, f_final, hist_final] = q3_jade(fun_fine, lb_f, ub_f, opts_f);
T_final = -f_final;
fprintf('========== 第4层结果：三导弹遮蔽时长之和 = %.4f s ==========\n', T_final);

%% ================== 6. 结果解析 ==================
nS_all = size(ginfo.slots, 1);
P_l = zeros(3, nS_all);  P_d = zeros(3, nS_all);
t_d = zeros(1, nS_all);  t_l = zeros(1, nS_all);
T_ind = zeros(1, nS_all);                    % 各弹位独立贡献时长
for s = 1:nS_all
    m_idx = ginfo.slots(s, 1);  uav_idx = ginfo.slots(s, 2);
    th = theta_v_all(uav_idx, 1);  v = theta_v_all(uav_idx, 2);
    u_h = [cosd(th); sind(th); 0];
    uav0 = params.uav0_all(:, uav_idx);
    t_l(s) = x_final(2*s - 1);  tau_s = x_final(2*s);
    P_l(:, s) = uav0 + v*t_l(s)*u_h;
    P_d(:, s) = P_l(:, s) + v*tau_s*u_h - [0; 0; 0.5*params.g*tau_s^2];
    t_d(s) = t_l(s) + tau_s;
end

% 各弹位独立贡献 + 各导弹独立遮蔽并集（主指标 = 三弹之和）
T_m = zeros(1, 3);  shielded_m = cell(1, 3);   % 各导弹遮蔽布尔序列
for m = 1:3
    rows = find(ginfo.slots(:, 1) == m);
    M0 = params.M0_all(:, m);  u_m = params.u_m_all(:, m);
    T_end = params.T_end(m);
    ts = 0 : params.dt : T_end;
    M_all = M0*ones(1, numel(ts)) + params.v_m*(u_m*ts);
    shielded = false(1, numel(ts));
    for r = rows(:)'
        mask = (ts >= t_d(r)) & (ts <= t_d(r) + 20);
        if any(mask)
            Cw = P_d(:, r)*ones(1, sum(mask)) - [0; 0; 3]*(ts(mask) - t_d(r));
            occ = q3_occlusion_window(M_all(:, mask), Cw, params.R, params);
            T_ind(r) = sum(occ)*params.dt;     % 该弹独立贡献
            shielded(mask) = shielded(mask) | occ(:)';
        end
    end
    T_m(m) = sum(shielded)*params.dt;
    shielded_m{m} = shielded;
end

T_total_final = sum(T_m);              % 主指标：三导弹遮蔽时长之和

%% ================== 7. 控制台输出 ==================
fprintf('\n========== 问题5 最优投放策略 ==========\n');
for s = 1:nS_all
    m_idx = ginfo.slots(s, 1);  uav_idx = ginfo.slots(s, 2);  slot = ginfo.slots(s, 3);
    fprintf('%s 弹%d -> %s：θ=%.2f°, v=%.2f m/s, t_l=%.3f s, tau=%.3f s, t_d=%.3f s\n', ...
        uav_names{uav_idx}, slot, m_names{m_idx}, ...
        theta_v_all(uav_idx, 1), theta_v_all(uav_idx, 2), t_l(s), x_final(2*s), t_d(s));
    fprintf('    投放点(%.1f, %.1f, %.1f)，起爆点(%.1f, %.1f, %.1f)，独立贡献 %.3f s\n', ...
        P_l(1, s), P_l(2, s), P_l(3, s), P_d(1, s), P_d(2, s), P_d(3, s), T_ind(s));
    if T_ind(s) < 0.1
        fprintf('    [提示] 该弹位与同机其他弹位窗口时段重叠，未能贡献额外遮蔽，结果表中留空。\n');
    end
end
fprintf('----------------------------------------------------------\n');
for m = 1:3
    fprintf('M%d 遮蔽时长（并集） = %.4f s\n', m, T_m(m));
end
fprintf('三导弹遮蔽时长之和（主指标） = %.4f s\n', T_total_final);
fprintf('==========================================================\n');

%% ================== 8. 写入 result3.xlsx ==================
header = {'无人机编号', '无人机运动方向', '无人机运动速度 (m/s)', ...
    '烟幕干扰弹编号', '烟幕干扰弹投放点的x坐标 (m)', ...
    '烟幕干扰弹投放点的y坐标 (m)', '烟幕干扰弹投放点的z 坐标 (m)', ...
    '烟幕干扰弹起爆点的x坐标 (m)', '烟幕干扰弹起爆点的y坐标 (m)', ...
    '烟幕干扰弹起爆点的z坐标 (m)', '有效干扰时长 (s)', '干扰的导弹编号'};
out = cell(18, 12);                     % 18 行，与模板 result3.xlsx 一致
out(1, :) = header;
for uav_idx = 1:5
    for slot = 1:3
        row = (uav_idx - 1)*3 + slot + 1;
        out(row, 1) = uav_names(uav_idx);
        out(row, 4) = {slot};
        s = find(ginfo.slots(:, 2) == uav_idx & ginfo.slots(:, 3) == slot, 1);
        % 仅填写有效弹位（独立贡献时长 >= 0.1s）；无效弹位留空：
        % 同一无人机同航向/速度下的多弹窗口时段由 θ/v 决定（物理预扫描
        % 结论），多余弹位与已生效弹位的窗口时段重叠，无法贡献额外遮蔽。
        if ~isempty(s) && T_ind(s) >= 0.1
            out(row, 2)  = {theta_v_all(uav_idx, 1)};
            out(row, 3)  = {theta_v_all(uav_idx, 2)};
            out(row, 5)  = {P_l(1, s)};
            out(row, 6)  = {P_l(2, s)};
            out(row, 7)  = {P_l(3, s)};
            out(row, 8)  = {P_d(1, s)};
            out(row, 9)  = {P_d(2, s)};
            out(row, 10) = {P_d(3, s)};
            out(row, 11) = {T_ind(s)};
            out(row, 12) = {m_names{ginfo.slots(s, 1)}};
        end
    end
end
out(18, 2) = {'注：以x轴为正向，逆时针方向为正，取值0~360（度）。'};
data_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
out_file = fullfile(data_dir, 'result3.xlsx');
writecell(out, out_file);
fprintf('结果已写入：%s\n', out_file);

%% ================== 9. 收敛曲线 ==================
figure('Color', 'w', 'Position', [80 80 1000 400]);
subplot(1, 3, 1);
hold on;
for m = 1:3
    if isempty(best_info.hist_lower{m}), continue; end
    plot(1:numel(best_info.hist_lower{m}), -best_info.hist_lower{m}, 'LineWidth', 1.3);
end
grid on; legend('M1', 'M2', 'M3', 'Location', 'best');
xlabel('进化代数'); ylabel('第2层各导弹遮蔽时长 (s)');
title('第2层 单导弹子问题收敛曲线');
subplot(1, 3, 2);
plot(1:numel(hist_final), -hist_final, 'r-', 'LineWidth', 1.5); grid on;
xlabel('进化代数'); ylabel('三导弹遮蔽时长之和 (s)');
title('第3/4层 联合精化收敛曲线');
subplot(1, 3, 3);
bar(T_m, 0.55);
set(gca, 'XTickLabel', {'M1', 'M2', 'M3'});
ylabel('遮蔽时长 (s)'); grid on;
title(sprintf('各导弹遮蔽时长（之和 = %.2f s）', T_total_final));

%% ================== 10. 遮蔽时间线图 ==================
figure('Color', 'w', 'Position', [100 100 1000 600]);
for m = 1:3
    subplot(3, 1, m);
    ts_m = 0 : params.dt : params.T_end(m);
    sh_m = false(1, numel(ts_m));
    rows = find(ginfo.slots(:, 1) == m);
    hold on; grid on;
    for r = rows(:)'
        mask = (ts_m >= t_d(r)) & (ts_m <= t_d(r) + 20);
        if any(mask)
            Cw = P_d(:, r)*ones(1, sum(mask)) - [0; 0; 3]*(ts_m(mask) - t_d(r));
            occ = q3_occlusion_window(params.M0_all(:, m)*ones(1, sum(mask)) + ...
                params.v_m*(params.u_m_all(:, m)*ts_m(mask)), Cw, params.R, params);
            sh_m(mask) = sh_m(mask) | occ(:)';
        end
    end
    stairs(ts_m, sh_m, 'k-', 'LineWidth', 1.5);
    for r = rows(:)'
        xline(t_d(r), '--', sprintf('%s弹%d', uav_names{ginfo.slots(r, 2)}, ginfo.slots(r, 3)));
    end
    xline(params.T_end(m), '--r', '到假目标');
    ylim([-0.1 1.1]);  yticks([0 1]);
    ylabel(sprintf('%s', m_names{m}));
    title(sprintf('%s 遮蔽时间线（并集 = %.2f s）', m_names{m}, T_m(m)));
    hold off;
end
xlabel('时间 t (s)');

%% ================== 11. 三维态势图 ==================
figure('Color', 'w', 'Position', [120 120 1000 760]);
hold on; grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('问题5 最优投放策略三维态势图');

colors_uav = {[0.2 0.6 0.2], [0 0.45 0.8], [0.85 0.55 0], [0.5 0.3 0.7], [0.9 0.2 0.3]};
colors_m = {'r', 'b', 'm'};

% 三枚导弹轨迹
for m = 1:3
    mt = linspace(0, params.T_end(m), 300);
    Mtraj = params.M0_all(:, m)' + params.v_m*mt'*params.u_m_all(:, m)';
    plot3(Mtraj(:, 1), Mtraj(:, 2), Mtraj(:, 3), colors_m{m}, 'LineWidth', 1.5);
    scatter3(params.M0_all(1, m), params.M0_all(2, m), params.M0_all(3, m), ...
             80, colors_m{m}, '^', 'filled');
    text(params.M0_all(1, m)+300, params.M0_all(2, m), params.M0_all(3, m), ...
         m_names{m}, 'Color', colors_m{m}, 'FontSize', 11);
end

% 各无人机：轨迹、投放点、起爆点、云团下沉
uav_used = unique(ginfo.slots(:, 2));
for u = 1:numel(uav_used)
    uav_idx = uav_used(u);
    u_h = [cosd(theta_v_all(uav_idx, 1)); sind(theta_v_all(uav_idx, 1)); 0];
    v_u = theta_v_all(uav_idx, 2);
    rows_u = find(ginfo.slots(:, 2) == uav_idx);
    t_max = max(t_l(rows_u));
    ut = linspace(0, t_max, 60);
    UAVtraj = params.uav0_all(:, uav_idx)' + v_u*ut'*u_h';
    plot3(UAVtraj(:, 1), UAVtraj(:, 2), UAVtraj(:, 3), '-', ...
          'Color', colors_uav{uav_idx}, 'LineWidth', 1.4);
    scatter3(params.uav0_all(1, uav_idx), params.uav0_all(2, uav_idx), ...
             params.uav0_all(3, uav_idx), 70, colors_uav{uav_idx}, 'o', 'filled');
    text(params.uav0_all(1, uav_idx)+300, params.uav0_all(2, uav_idx), ...
         params.uav0_all(3, uav_idx), uav_names{uav_idx}, ...
         'Color', colors_uav{uav_idx}, 'FontSize', 11);
    for r = rows_u(:)'
        scatter3(P_l(1, r), P_l(2, r), P_l(3, r), 80, 'g', 'd', 'filled');
        scatter3(P_d(1, r), P_d(2, r), P_d(3, r), 100, 'k', 'p', 'filled');
        cdt = linspace(t_d(r), min(t_d(r)+20, params.T_end(ginfo.slots(r, 1))), 50);
        Ctraj = P_d(:, r) - [0; 0; 3]*(cdt - t_d(r));
        plot3(Ctraj(1, :), Ctraj(2, :), Ctraj(3, :), 'Color', [0.6 0.6 0.6], ...
              'LineWidth', 1.0);
        text(P_l(1, r)+200, P_l(2, r), P_l(3, r), ...
             sprintf('%s弹%d', uav_names{uav_idx}, ginfo.slots(r, 3)), ...
             'Color', 'g', 'FontSize', 8);
    end
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
view(3);  hold off;

disp('求解完成。请检查控制台输出与 result3.xlsx，以及三张图。');

%% ================== 局部函数：第2层种子生成 ==================
function seeds = q5_build_seeds(m_idx, servers, scan_par, params, par_relay)
% 为单导弹子问题（第2层）构造初始种群。
% 决策变量 = [θ/v (2*nU), t_l/τ (2*nS)]，前 2*nU 维为各服务机 [θ,v]。
% 种子设计：
%   1) 前 10 组：预扫描单弹 θ/v + 接力速度档，t_l 沿弹位错开 1 s；
%   2) 追加接力预扫描最优解（par_relay），使 JADE 从接力盆地出发，
%      修正「单弹最优 v/τ 不产生接力增益、JADE 因而漂移」的问题。
nU = numel(servers.uavs);
nS = size(servers.slots, 1);
D  = 2*nU + 2*nS;

% 每服务机的预扫描单弹 [θ, v, t_l, τ]
base = zeros(nU, 4);
for u = 1:nU
    uav_idx = servers.uavs(u);
    r = find(scan_par(:, 1) == uav_idx & scan_par(:, 2) == m_idx, 1);
    if isempty(r)
        base(u, :) = [180, 110, 5, 4];      % 兜底
    else
        base(u, :) = scan_par(r, 3:6);      % [θ, v, t_l, τ]
    end
end

% 服务机编号 -> 在 servers.uavs 中的位置
uav_pos = zeros(1, max(servers.uavs));
for u = 1:nU
    uav_pos(servers.uavs(u)) = u;
end

% 接力速度档（数值验证：接力方向 = 单弹最优方向，多弹沿该方向 1 s 错开）
v_relay = [85, 90, 95, 100, 105];

seeds = zeros(10, D);
for g = 1:10
    % ---- θ/v 部分 ----
    for u = 1:nU
        th = base(u, 1);                     % 单弹最优 θ 即接力方向
        if g <= 5
            v = base(u, 2);                  % 预扫描单弹最优 v
        else
            v = v_relay(g-5);                % 不同速度档探索接力
        end
        seeds(g, 2*u-1) = th;
        seeds(g, 2*u)   = v;
    end
    % ---- t_l/τ 部分：同机多弹沿航向错开 1 s（相邻投放间隔 >= 1 s）----
    for s = 1:nS
        uav_idx = servers.slots(s, 1);
        u  = uav_pos(uav_idx);
        q  = servers.slots(s, 2);            % 弹位编号 1~3
        tl  = base(u, 3) + (q-1)*1 + 0.1*(g-1);
        tau = base(u, 4) + 0.1*mod(g-1, 5);
        if tl < 0, tl = 0; end
        seeds(g, 2*nU + 2*s-1) = tl;
        seeds(g, 2*nU + 2*s)   = tau;
    end
end

% ---- 追加接力预扫描最优解作为种子 ----
if nargin >= 5 && ~isempty(par_relay)
    extra = zeros(0, D);
    for u = 1:nU
        uav_idx = servers.uavs(u);
        xr = par_relay{uav_idx, m_idx};
        if isempty(xr), continue; end
        row = zeros(1, D);
        % θ/v：本机用接力最优，其余机用单弹最优
        for uu = 1:nU
            row(2*uu-1) = base(uu, 1);
            row(2*uu)   = base(uu, 2);
        end
        row(2*u-1) = xr(1);
        row(2*u)   = xr(2);
        % t_l/τ：本机 3 弹用接力最优，其余机单弹最优错开 1 s
        for s = 1:nS
            si = servers.slots(s, 1);
            q  = servers.slots(s, 2);
            uu = uav_pos(si);
            if si == uav_idx
                row(2*nU + 2*s-1) = xr(2*q + 1);   % t_l
                row(2*nU + 2*s)   = xr(2*q + 2);   % τ
            else
                row(2*nU + 2*s-1) = base(uu, 3) + (q-1)*1;
                row(2*nU + 2*s)   = base(uu, 4);
            end
        end
        extra(end+1, :) = row; %#ok<AGROW>
    end
    seeds = [seeds; extra];
end
end
