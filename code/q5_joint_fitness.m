function f = q5_joint_fitness(x, params, ginfo, theta_v_all)
% q5_joint_fitness  问题5【第3层】三导弹联合精化目标（求和/并集主指标）
%
% 在固定每架无人机航向/速度（第2层结果，theta_v_all）的前提下，
% 对所有弹位的投放/延时 (t_l, tau) 做联合优化：
%  决策变量 x = [t_l1, tau1, ..., t_l_nS, tau_nS]   （2*nS 维）
% 目标：最大化三枚导弹【各自独立遮蔽并集时长之和】
%   f = -(T_M1 + T_M2 + T_M3) + 约束惩罚（最小化）
%
% 说明：
%   1) 主指标 = 三枚导弹各自被其指派弹位云团遮蔽的并集时长之和
%      （f_total = f_M1 + f_M2 + f_M3），与问题5求解思路（second）
%      一致；每枚导弹在其完整飞行时段 [0, T_end(m)] 上独立统计。
%   2) 联合优化的意义：全部弹位统一搜索，显式处理同一无人机多弹位
%      的连投间隔约束，避免子问题分解后跨导弹资源冲突。
%
% 约束（惩罚，w=1e4）：
%   1) 各弹起爆 t_d <= 对应导弹到达假目标时刻；
%   2) 各弹 tau <= 落地前起爆上限；
%   3) 同一无人机相邻两弹投放间隔 >= 1 s。
%
% 输入：
%   x           - 决策变量（2*nS 维）
%   params      - 场景结构体
%   ginfo       - 全局弹位表：.slots = nS x 3（每行 [导弹编号, 无人机编号,
%                弹位编号]），行序与 x 一一对应
%   theta_v_all - 5 x 2 每架无人机固定的航向角/速度（第2层结果）
%
% 依赖：q3_occlusion_window.m

nS = size(ginfo.slots, 1);
w  = 1e4;  pen = 0;

% ---- 各导弹遮蔽并集时长（各自完整时间轴）----
T_total = 0;
for m = 1:3
    rows = find(ginfo.slots(:, 1) == m);
    if isempty(rows), continue; end

    M0    = params.M0_all(:, m);
    u_m   = params.u_m_all(:, m);
    T_end = norm(M0)/params.v_m;
    ts = 0 : params.dt : T_end;
    nT = numel(ts);
    M_all = M0*ones(1, nT) + params.v_m*(u_m*ts);
    shielded = false(1, nT);

    for r = rows(:)'
        t_l = x(2*r-1);  tau = x(2*r);
        uav_idx = ginfo.slots(r, 2);
        th = theta_v_all(uav_idx, 1);  v = theta_v_all(uav_idx, 2);
        u_h = [cosd(th); sind(th); 0];
        uav0 = params.uav0_all(:, uav_idx);
        tau_max = sqrt(2*uav0(3)/params.g);

        P_l = uav0 + v*t_l*u_h;
        P_d = P_l + v*tau*u_h - [0; 0; 0.5*params.g*tau^2];
        t_d = t_l + tau;

        mask = (ts >= t_d) & (ts <= t_d + 20);
        if any(mask)
            Cw = P_d*ones(1, sum(mask)) - [0; 0; 3]*(ts(mask) - t_d);
            occ = q3_occlusion_window(M_all(:, mask), Cw, params.R, params);
            shielded(mask) = shielded(mask) | occ(:)';
        end

        if t_d > T_end,       pen = pen + w*(t_d - T_end)^2;   end
        if tau > tau_max,     pen = pen + w*(tau - tau_max)^2; end
    end
    T_total = T_total + sum(shielded)*params.dt;
end

% ---- 同机连投间隔 >= 1 s ----
for a = 1:nS
    for b = a+1 : nS
        if ginfo.slots(a, 2) == ginfo.slots(b, 2)
            dgap = abs(x(2*a-1) - x(2*b-1));
            if dgap < 1
                pen = pen + w*(1 - dgap)^2;
            end
        end
    end
end

f = -T_total + pen;
end
