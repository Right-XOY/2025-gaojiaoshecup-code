function f = q4_single_fitness(x, params, uav0, target_center)
% q4_single_fitness  问题4【下层】单机投放参数优化目标函数（4 维软约束）
%
% 在问题4的分层优化中，下层固定上层给定的时序目标，对单架无人机独立优化：
%   决策变量 x = [theta_deg, v, t_l, tau]
%     theta_deg - 该机航向角（度，0~360，x 轴正向逆时针为正）
%     v         - 该机飞行速度 (m/s)，70~140
%     t_l       - 烟幕弹投放时刻 (s)
%     tau       - 引信延时 (s)，落地前起爆
%
% 目标：软约束形式的"最大化自身遮蔽时长 + 窗口中心贴近目标时刻"
%   f = -T_occ + w_center*(窗口中心 - target_center)^2 + 约束惩罚
%   - T_occ        : 该枚烟幕弹对 M1 的独立有效遮蔽时长（并集只含单弹）
%   - 窗口中心     : 该云团遮蔽时间段的起点终点均值，软约束允许遮蔽
%                    窗口"延伸"出上层分配区间之外，只追求中心对齐
%   - 约束惩罚     : 起爆须在导弹到达假目标前；引信延时须在落地前
%
% 物理过程（与问题2/3完全一致）：
%   投放点 P_l = uav0 + v*t_l*u_h
%   起爆点 P_d = P_l + v*tau*u_h - [0,0,0.5*g*tau^2]
%   云团   C(t)= P_d - [0,0,3]*(t - t_d)，t in [t_d, t_d+20]
%
% 依赖：q3_occlusion_window.m

theta_deg = x(1);  v = x(2);  t_l = x(3);  tau = x(4);
u_h = [cosd(theta_deg); sind(theta_deg); 0];          % 水平航向单位向量
h   = uav0(3);                                        % 无人机飞行高度 (m)
tau_max = sqrt(2*h/params.g);                         % 落地前起爆的引信延时上限

% ---- 时间网格与导弹轨迹 ----
ts = 0 : params.dt : params.T_end;
M_all = params.M0*ones(1, numel(ts)) + params.v_m*(params.u_m*ts);

% ---- 单枚弹：投放点、起爆点、云团窗口 ----
P_l = uav0 + v*t_l*u_h;
P_d = P_l + v*tau*u_h - [0; 0; 0.5*params.g*tau^2];
t_d = t_l + tau;                                      % 起爆时刻

T_occ = 0;  center_pen = 0;
idx = find(ts >= t_d & ts <= t_d + 20);
if ~isempty(idx)
    Cw = P_d*ones(1, numel(idx)) - [0; 0; 3]*(ts(idx) - t_d);
    occ = q3_occlusion_window(M_all(:, idx), Cw, params.R, params);
    T_occ = sum(occ)*params.dt;
    if any(occ)
        k1 = find(occ, 1, 'first');  k2 = find(occ, 1, 'last');
        center = (ts(idx(k1)) + ts(idx(k2)))/2;       % 遮蔽窗口中心
        center_pen = params.w_center*(center - target_center)^2;
    else
        % 无遮蔽时以云团有效期中点 (t_d+10) 作参考，引导搜索靠近目标时刻
        center_pen = params.w_center*(t_d + 10 - target_center)^2;
    end
end

% ---- 约束惩罚 ----
w = 1e4;  pen = 0;
if t_d > params.T_end                                % 起爆须在导弹到达假目标前
    pen = pen + w*(t_d - params.T_end)^2;
end
if tau > tau_max                                     % 引信延时落地前起爆
    pen = pen + w*(tau - tau_max)^2;
end

f = -T_occ + center_pen + pen;
end
