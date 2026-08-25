function f = q3_fitness(x, params)
% q3_fitness  问题3目标函数：FY1 投放 3 枚烟幕干扰弹干扰 M1，最大化总遮蔽时长
%
% 决策变量 x = [theta_deg, v, t_l1, tau1, t_l2, tau2, t_l3, tau3]
%   theta_deg - FY1 航向角（度，0~360，x 轴正向逆时针为正）
%   v         - FY1 飞行速度 (m/s)，70~140
%   t_li      - 第 i 枚弹投放时刻 (s)，要求 t_l2-t_l1>=1、t_l3-t_l2>=1
%   taui      - 第 i 枚弹引信延时 (s)，落地前起爆
%
% 物理过程：
%   投放点  P_l,i = P0 + v*t_li*u_h
%   起爆点  P_d,i = P_l,i + v*taui*u_h - [0,0,0.5*g*taui^2]
%   云团 i  C_i(t) = P_d,i - [0,0,3*(t-t_di)]，t in [t_di, t_di+20]
%   导弹    M(t) = M0 + 300*t*u_m
%
% 目标值：总有效遮蔽时长 = 三个云团遮蔽区间【并集】的长度（s），
%        任一时刻任一枚云团完全遮蔽圆柱即计为遮蔽。
%
% 约束（惩罚函数处理）：
%   1) 投放间隔：t_l2-t_l1>=1, t_l3-t_l2>=1
%   2) 起爆须在导弹到达假目标前：t_di = t_li+taui <= T_end
% 返回负遮蔽时长（供最小化）。

theta_deg = x(1);  v = x(2);
t_l = x(3:2:7);    % [t_l1 t_l2 t_l3]
tau = x(4:2:8);    % [tau1 tau2 tau3]

u_h = [cosd(theta_deg); sind(theta_deg); 0];     % 无人机水平航向单位向量

% ---- 时间网格与导弹轨迹（一次计算，三枚弹共用）----
ts = 0 : params.dt : params.T_end;
T  = numel(ts);
M_all = params.M0*ones(1, T) + params.v_m*(params.u_m*ts);   % 3 x T

% ---- 并集遮蔽时长 ----
shielded = false(1, T);
for i = 1:3
    P_l = params.uav0 + v*t_l(i)*u_h;                        % 投放点
    P_d = P_l + v*tau(i)*u_h - [0; 0; 0.5*params.g*tau(i)^2];% 起爆点
    t_d = t_l(i) + tau(i);                                   % 起爆时刻
    mask = (ts >= t_d) & (ts <= t_d + 20);                   % 云团有效窗口（逻辑掩码）
    if any(mask)
        Cw = P_d*ones(1, sum(mask)) - [0; 0; 3]*(ts(mask) - t_d);
        occ = q3_occlusion_window(M_all(:, mask), Cw, params.R, params);
        if numel(occ) ~= sum(mask)                           % 防御：形状不匹配立刻指明
            error('q3_occlusion_window 返回 %d 个元素，期望 %d 个——请 clear functions 重载后重试', ...
                  numel(occ), sum(mask));
        end
        shielded(mask) = shielded(mask) | occ(:)';           % 两侧均行向量（逻辑索引保留行方向）
    end
end
T_total = sum(shielded)*params.dt;

% ---- 约束惩罚 ----
w = 1e4;
pen = 0;
for i = 1:2
    gap = t_l(i+1) - t_l(i) - 1;               % 投放间隔须 >= 1s
    if gap < 0, pen = pen + w*gap^2; end
end
for i = 1:3
    if t_l(i) + tau(i) > params.T_end          % 起爆须在导弹到达假目标前
        pen = pen + w*(t_l(i) + tau(i) - params.T_end)^2;
    end
end

f = -T_total + pen;
end
