function f = q2_fitness(x, params)
% q2_fitness  目标函数：给定投放策略，计算对 M1 的有效遮蔽时长（取负值供最小化）
%
% 决策变量 x = [theta, v, t_l, tau]
%   theta - FY1 航向角（水平面，rad），0~2*pi
%   v     - FY1 飞行速度 (m/s)，70~140
%   t_l   - 受领任务后到投放干扰弹的时间 (s)
%   tau   - 引信延时，投放后到起爆的时间 (s)
%
% 输出：
%   f     - 负的有效遮蔽时长 (s)，供 ga/patternsearch 最小化
%
% 物理过程：
%   投放点  P_l = P0 + v*t_l*(cos,sin,0)
%   起爆点  P_d = P_l + v*tau*(cos,sin,0) - [0,0,0.5*g*tau^2]
%   云团    C(t) = P_d - [0,0,3*(t-t_det)]，t in [t_det, t_det+20]
%   导弹    M(t) = M0 + 300*t*u_m
%   遮蔽时长：在 [t_det, t_det+20] ∩ [0, T_end] 内，
%            圆柱目标被完全遮蔽的时间长度（T_end 为导弹到达假目标时刻）

theta = x(1);
v     = x(2);
t_l   = x(3);
tau   = x(4);
g     = params.g;

u_h = [cos(theta); sin(theta); 0];          % 无人机水平航向单位向量

P_l = params.uav0 + v*t_l*u_h;              % 投放点
P_d = P_l + v*tau*u_h - [0; 0; 0.5*g*tau^2];% 起爆点
t_det = t_l + tau;                          % 起爆时刻

T_shield = 0;
% 起爆须在导弹到达假目标之前，且起爆点不低于地面
if t_det < params.T_end && P_d(3) >= 0
    dt     = params.dt;
    t_stop = min(t_det + 20, params.T_end); % 云团有效窗口与交战窗口取交集
    t0     = t_det + dt/2;                  % 中点采样，减小积分误差
    for t = t0 : dt : t_stop
        C = P_d - [0; 0; 3*(t - t_det)];    % 云团匀速下沉 3 m/s
        M = params.M0 + 300*t*params.u_m;   % 导弹位置
        if q2_occlusion(M, C, params.R, params)
            T_shield = T_shield + dt;
        end
    end
end

f = -T_shield;                              % 最小化负遮蔽时长 = 最大化遮蔽时长
end
