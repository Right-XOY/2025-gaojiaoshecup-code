function [U, delta] = q5_pair_prescan(params, pot_relay, par_relay)
% q5_pair_prescan  问题5【两机协同预扫描】对每个可行的两机组合 (j1<j2, m)，
%   用 JADE 求「两机协同服务 m」的并集最优遮蔽时长 U(j1,j2,m)，
%   并计算净重叠损失 delta = U - pot_relay(j1,m) - pot_relay(j2,m) (<=0)。
%
% 目的：修正第1层 ILP「相加」目标高估多机协同的问题。两机同服务一枚
%   导弹时遮蔽窗口会重叠，真实并集 <= 两机单独接力时长之和；用 delta
%   修正 ILP 目标，使候选排名逼近真实并集（而非简单相加）。
%
% 实现：复用 q5_missile_fitness（两机 6 弹，决策变量 16 维 =
%   [θ1,v1,θ2,v2, 6 弹各自 t_l/τ]），以两机各自的单机接力最优解 par_relay
%   拼接为种子，JADE 精化求两机协同并集。
%
% 输入：params（场景，含 dt）、pot_relay（单机接力 5x3）、
%       par_relay（单机接力最优解 cell，每项 8 维 [θ,v,tl1,τ1,tl2,τ2,tl3,τ3]）
% 输出：U     - 5x5x3 数组，U(j1,j2,m) 两机并集最优时长（j1<j2 处有值，其余 0）
%       delta - 5x5x3 数组，delta(j1,j2,m) = U - pot_relay(j1,m) - pot_relay(j2,m)

nU = 5; nM = 3;
U = zeros(nU, nU, nM);
delta = zeros(nU, nU, nM);

opts.NP = 40; opts.MAXITER = 200; opts.p = 0.05; opts.c = 0.1;

for m = 1:nM
    feas = find(pot_relay(:, m) > 0)';        % 该导弹的可行服务机
    nf = numel(feas);
    for ii = 1:nf
        for jj = ii+1:nf
            j1 = feas(ii); j2 = feas(jj);

            % ---- 两机 6 弹：机1 弹1~3 在前，机2 弹1~3 在后 ----
            servers.uavs  = [j1, j2];
            servers.slots = [j1 1; j1 2; j1 3; j2 1; j2 2; j2 3];

            % ---- 边界：[θ1,v1,θ2,v2] + 6 弹 [t_l, τ] ----
            h1 = params.uav0_all(3, j1); tau_max1 = sqrt(2*h1/params.g);
            h2 = params.uav0_all(3, j2); tau_max2 = sqrt(2*h2/params.g);
            T_end = params.T_end(m);
            lb = [0, params.v_min, 0, params.v_min, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0];
            ub = [360, params.v_max, 360, params.v_max, ...
                  T_end, tau_max1, T_end, tau_max1, T_end, tau_max1, ...
                  T_end, tau_max2, T_end, tau_max2, T_end, tau_max2];

            % ---- 种子：两机各自的单机接力最优解拼接 ----
            p1 = par_relay{j1, m};  p2 = par_relay{j2, m};
            x0 = zeros(1, 16);
            x0(1:2)  = p1(1:2);                 % θ1, v1
            x0(3:4)  = p2(1:2);                 % θ2, v2
            x0(5:2:9)   = p1(3:2:7);            % 机1 tl1,tl2,tl3
            x0(6:2:10)  = p1(4:2:8);            % 机1 τ1,τ2,τ3
            x0(11:2:15) = p2(3:2:7);            % 机2 tl1,tl2,tl3
            x0(12:2:16) = p2(4:2:8);            % 机2 τ1,τ2,τ3

            % 变体：整体平移某机投放时刻，探索两机窗口错开
            S = x0;
            for d = [0.5 1 2]
                xv = x0; xv(11:2:15) = xv(11:2:15) + d;   % 机2 整体延后
                S(end+1, :) = xv; %#ok<AGROW>
            end
            for d = [0.5 1]
                xv = x0; xv(5:2:9) = xv(5:2:9) + d;       % 机1 整体延后
                S(end+1, :) = xv; %#ok<AGROW>
            end
            opts.X0 = S;

            fun = @(x) q5_missile_fitness(x, params, m, servers);
            [~, fbest] = q3_jade(fun, lb, ub, opts);

            U(j1, j2, m) = max(-fbest, max(pot_relay(j1, m), pot_relay(j2, m)));
            delta(j1, j2, m) = U(j1, j2, m) - pot_relay(j1, m) - pot_relay(j2, m);
        end
    end
end

fprintf('两机协同预扫描完成（delta = 两机并集 - 两机单机接力之和，<=0）。\n');
end
