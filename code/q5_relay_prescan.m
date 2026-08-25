function [pot_relay, par_relay] = q5_relay_prescan(params, pot, scan_par)
% q5_relay_prescan  问题5【接力预扫描】对每个可行 (无人机, 导弹) 组合，
%   用 JADE 求「单机 3 弹接力」的最优遮蔽时长（θ/v 共享 + 每弹独立 t_l/τ）。
%
% 目的（修正旧假设「同机多弹无法接力」）：
%   1) 接力时长 pot_relay 作为第1层 ILP 系数。单机 3 弹沿单弹最优方向
%      （迎向导弹）接力可显著延长遮蔽（FY2→M2 3.9→7.5s、FY1→M1 4.6→6.0s），
%      据此指派能避免把多余机（如 FY4）错误地压到已被单机接力覆盖的
%      导弹上造成浪费。
%   2) par_relay 作为第2层 JADE 的接力种子，使接力盆地被可靠发现
%      （旧种子只用单弹最优 v/τ，该配置不产生接力增益，JADE 因而漂移）。
%
% 关键物理结论（数值验证）：
%   接力方向 = 单弹最优方向（迎向导弹），而非其反方向；
%   接力需特定 v 与每弹各自 τ，单弹最优 v/τ 未必能接力。
%
% 输入：params（场景，需已设 params.dt）、pot（单弹潜力 5x3）、
%       scan_par（单弹最优参数表，[uav,m,θ,v,t_l,τ,...]）
% 输出：pot_relay（接力时长 5x3）、par_relay（cell，每项 8 维最优解）

nU = 5; nM = 3;
pot_relay = zeros(nU, nM);
par_relay = cell(nU, nM);

opts.NP = 40;  opts.MAXITER = 200;  opts.p = 0.05;  opts.c = 0.1;

for uav_idx = 1:nU
    h = params.uav0_all(3, uav_idx);
    tau_max = sqrt(2*h/params.g);
    for m = 1:nM
        if pot(uav_idx, m) <= 0, continue; end
        r = find(scan_par(:,1)==uav_idx & scan_par(:,2)==m, 1);
        if isempty(r), continue; end
        th0  = scan_par(r,3);
        v0   = scan_par(r,4);
        tl0  = scan_par(r,5);
        tau0 = scan_par(r,6);

        servers.uavs  = uav_idx;
        servers.slots = [uav_idx 1; uav_idx 2; uav_idx 3];

        T_end = params.T_end(m);
        lb = [0, params.v_min, 0,0, 0,0, 0,0];
        ub = [360, params.v_max, T_end, tau_max, T_end, tau_max, T_end, tau_max];

        % ---- 多样种子：覆盖 v、τ（同 τ 网格）----
        v_vals   = params.v_min:10:params.v_max;
        tau_vals = unique([0.3 0.5 1 2 4 tau0]);
        tau_vals = tau_vals(tau_vals > 0 & tau_vals <= tau_max);
        S = zeros(0, 8);
        for vv = v_vals
            for tt = tau_vals
                tls = max([tl0, tl0+1, tl0+2], 0);
                if tls(3) + tt > T_end, continue; end
                S(end+1,:) = [th0, vv, tls(1), tt, tls(2), tt, tls(3), tt]; %#ok<AGROW>
            end
        end
        % ---- 不同 τ 接力变体（弹2 用更短 τ，使云团落在更晚/更近导弹路径段）----
        for dtau = [0.5 1 2]
            tt1 = tau0;  tt2 = max(0.3, tau0 - dtau);
            tls = max([tl0, tl0+1, tl0+2], 0);
            if tls(3) + max(tt1,tt2) > T_end, continue; end
            S(end+1,:) = [th0, v0, tls(1), tt1, tls(2), tt2, tls(3), tt1]; %#ok<AGROW>
        end

        if isempty(S)
            pot_relay(uav_idx,m) = pot(uav_idx,m);
            continue;
        end
        if size(S,1) > opts.NP, S = S(1:opts.NP,:); end
        opts.X0 = S;

        fun = @(x) q5_missile_fitness(x, params, m, servers);
        [xbest, fbest] = q3_jade(fun, lb, ub, opts);
        pot_relay(uav_idx, m) = max(-fbest, pot(uav_idx, m));
        par_relay{uav_idx, m} = xbest;
    end
end

fprintf('接力预扫描完成：pot_relay 矩阵（行=FY1..FY5，列=M1..M3）：\n');
disp(pot_relay);
end
