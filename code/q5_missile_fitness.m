function f = q5_missile_fitness(x, params, m_idx, servers)
% q5_missile_fitness  问题5【第2层】单枚导弹多弹接力遮蔽子问题目标
%
% 决策变量 x = [θ_1,v_1, ..., θ_nU,v_nU,  t_l1,τ_1, ..., t_l_nS,τ_nS]
%   前 2*nU 维：各服务机的航向角/速度（同机多弹共享，满足"θ/v 一旦
%               确定不再调整"的题目约束；θ∈[0,360], v∈[70,140]）
%   后 2*nS 维：各弹位投放时刻/引信延时
% 目标：最大化该导弹被所有云团遮蔽的【独立并集时长】
%   f = -T_occ + 约束惩罚（最小化）
%
% 说明：θ/v 纳入决策变量后，JADE 可在"单弹最优方向（朝+x）"与
%   "接力方向（朝-x 迎向导弹）"之间自行权衡，使同机多弹真正形成接力。
%
% 约束（惩罚，w=1e4）：
%   1) 各弹起爆 t_d = t_l + tau <= 该导弹到达假目标时刻 T_end_m；
%   2) 各弹引信延时 tau <= 落地前起爆上限 sqrt(2*h/g)；
%   3) 同一架无人机的相邻两弹投放时刻至少间隔 1 s。
%
% 输入：
%   x       - 决策变量（2*nU + 2*nS 维）
%   params  - 场景结构体（含 M0_all 3x3、u_m_all 3x3、uav0_all 3x5 等）
%   m_idx   - 导弹编号 1~3
%   servers - 结构体：
%             .uavs  = 1 x nU 服务机编号（去重）
%             .slots = nS x 2 弹位表（每行 [无人机编号, 弹位编号]）
%
% 依赖：q3_occlusion_window.m

M0   = params.M0_all(:, m_idx);            % 导弹初始位置
u_m  = params.u_m_all(:, m_idx);           % 导弹飞行方向
T_end= norm(M0)/params.v_m;                % 到达假目标时刻

ts = 0 : params.dt : T_end;
M_all = M0*ones(1, numel(ts)) + params.v_m*(u_m*ts);

nU = numel(servers.uavs);
nS = size(servers.slots, 1);

% ---- 解析决策变量 ----
thv = x(1 : 2*nU);                         % 各服务机 [θ, v] 交错排列
t_l = x(2*nU+1 : 2 : 2*nU+2*nS);           % 各弹位投放时刻
tau = x(2*nU+2 : 2 : 2*nU+2*nS);           % 各弹位引信延时

% 服务机编号 -> 其在 servers.uavs 中的位置（用于取 θ/v）
uav_pos = zeros(1, max(servers.uavs));
for u = 1:nU
    uav_pos(servers.uavs(u)) = u;
end

shielded = false(1, numel(ts));
w = 1e4;  pen = 0;

% ---- 各弹位遮蔽叠加（并集） ----
for s = 1:nS
    uav_idx = servers.slots(s, 1);
    u = uav_pos(uav_idx);
    th = thv(2*u-1);  v = thv(2*u);
    u_h = [cosd(th); sind(th); 0];
    uav0 = params.uav0_all(:, uav_idx);
    tau_max = sqrt(2*uav0(3)/params.g);

    P_l = uav0 + v*t_l(s)*u_h;                                % 投放点
    P_d = P_l + v*tau(s)*u_h - [0; 0; 0.5*params.g*tau(s)^2]; % 起爆点
    t_d = t_l(s) + tau(s);

    mask = (ts >= t_d) & (ts <= t_d + 20);                     % 云团有效期（逻辑掩码）
    if any(mask)
        Cw = P_d*ones(1, sum(mask)) - [0; 0; 3]*(ts(mask) - t_d);
        occ = q3_occlusion_window(M_all(:, mask), Cw, params.R, params);
        shielded(mask) = shielded(mask) | occ(:)';
    end

    if t_d > T_end,        pen = pen + w*(t_d - T_end)^2;     end
    if tau(s) > tau_max,   pen = pen + w*(tau(s) - tau_max)^2; end
end

% ---- 同机连投间隔 >= 1 s ----
for a = 1:nS
    for b = a+1 : nS
        if servers.slots(a, 1) == servers.slots(b, 1)
            dgap = abs(t_l(a) - t_l(b));
            if dgap < 1
                pen = pen + w*(1 - dgap)^2;
            end
        end
    end
end

T_occ = sum(shielded)*params.dt;
f = -T_occ + pen;
end
