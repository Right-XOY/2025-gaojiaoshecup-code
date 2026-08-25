function f = q4_joint_fitness(x, params)
% q4_joint_fitness  问题4【精化层】12 维联合并集目标函数
%
% 决策变量 x（按无人机分组，每组 4 维）：
%   x = [theta1,v1,tl1,tau1,  theta2,v2,tl2,tau2,  theta3,v3,tl3,tau3]
% 对应 FY1、FY2、FY3 各 1 枚烟幕弹。
%
% 目标：最大化三架无人机云团对 M1 遮蔽的【总并集时长】，
%       同一时刻只要任意一枚云团完全遮蔽圆柱即计入遮蔽。
%   f = -T_union + 约束惩罚
%
% 约束（惩罚函数处理）：
%   1) 各弹起爆须在导弹到达假目标前：t_di = t_li + taui <= T_end
%   2) 各弹引信延时须在落地前：taui <= sqrt(2*h_i/g)
%   3) 速度、航向边界由 JADE 的 lb/ub 保证
%
% 依赖：q3_occlusion_window.m

ts = 0 : params.dt : params.T_end;
T  = numel(ts);
M_all = params.M0*ones(1, T) + params.v_m*(params.u_m*ts);

shielded = false(1, T);
w = 1e4;  pen = 0;

for i = 1:3
    theta_deg = x(4*i-3);  v = x(4*i-2);  t_l = x(4*i-1);  tau = x(4*i);
    u_h = [cosd(theta_deg); sind(theta_deg); 0];
    uav0 = params.uav0_all(:, i);                          % 第 i 架无人机初始位置
    tau_max = sqrt(2*uav0(3)/params.g);                    % 该机落地前起爆上限

    P_l = uav0 + v*t_l*u_h;                                % 投放点
    P_d = P_l + v*tau*u_h - [0; 0; 0.5*params.g*tau^2];    % 起爆点
    t_d = t_l + tau;                                       % 起爆时刻

    mask = (ts >= t_d) & (ts <= t_d + 20);                  % 云团有效窗口（逻辑掩码）
    if any(mask)
        Cw = P_d*ones(1, sum(mask)) - [0; 0; 3]*(ts(mask) - t_d);
        occ = q3_occlusion_window(M_all(:, mask), Cw, params.R, params);
        shielded(mask) = shielded(mask) | occ(:)';
    end

    if t_l + tau > params.T_end, pen = pen + w*(t_l + tau - params.T_end)^2; end
    if tau > tau_max,            pen = pen + w*(tau - tau_max)^2; end
end

T_union = sum(shielded)*params.dt;
f = -T_union + pen;
end
