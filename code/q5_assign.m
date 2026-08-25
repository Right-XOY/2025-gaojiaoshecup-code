function [cand, nGet] = q5_assign(pot, k)
% q5_assign  问题5【第1层】任务分配——精确 pot + ILP 自动求解（不硬编码）
%
% 分配模型（ILP 框架，其中 c_ij 采用精确物理预扫描得到的 pot 矩阵，
% 即实际可达遮蔽时长，而非理想化上界 20 - T_fly）：
%     max  Σ_j Σ_m  pot(j,m) * y(j,m)
%     s.t. Σ_m y(j,m) <= 1        （每架无人机至多服务 1 种导弹：
%                                   一架机一个航向/速度只能遮蔽一枚
%                                   导弹的飞行路径）
%          Σ_j y(j,m) >= 1        （每枚导弹至少被 1 架无人机干扰）
%          y(j,m) ∈ {0,1}
%  说明：
%   1) 同机多弹可沿单弹最优方向形成接力（精确预扫描结论），故 pot 取
%      单机 3 弹接力最优时长（pot_relay），y(j,m) 取 0/1 表示该机是否
%      服务该导弹，弹位按 3 弹展开由第2层决定；
%   2) 每架无人机实际允许投放至多 3 枚烟幕弹，弹位展开由第2层优化
%      决定（同机多弹若窗口重叠会自动留空），候选表中弹数记为 3。
%   3) 求解方式：小规模整数规划直接枚举（5 机各 4 种选择，共 4^5
%      = 1024 组合），不依赖整数规划工具箱，也不硬编码任何分配。
%
% 输入：
%   pot - 5 x 3 遮蔽潜力矩阵（行=无人机 FY1~FY5，列=导弹 M1~M3），
%         单机 3 弹接力最优时长，由 q5_relay_prescan.m 生成（pot_relay）
%   k   - 候选上限（默认 3；有效候选不足时返回实际数目）
% 输出：
%   cand - 1 x nGet 结构体数组；cand(t).missile{m} = n x 2 服务表
%          （第1列无人机编号，第2列该机为此导弹分配的弹位上限）
%   nGet - 实际返回的候选套数 = min(k, 有效候选数)
%
% 依赖：无（纯指派层）

if nargin < 2, k = 3; end
nU = size(pot, 1);  nM = size(pot, 2);
base = nM + 1;                        % 每机选择：0=不出动, 1..nM=服务导弹

% ---- 枚举全部 4^nU 种分配 ----
codes = dec2base(0:base^nU-1, base) - '0';    % (base^nU) x nU，取值 0..nM
nC = size(codes, 1);
scores = zeros(nC, 1);
assigns = cell(nC, 1);
valid  = false(nC, 1);

for c = 1:nC
    a = codes(c, :);                  % 1 x nU：每机的服务对象
    covered = false(1, nM);
    sc = 0;  ok = true;
    for j = 1:nU
        m = a(j);
        if m == 0, continue; end      % 该机不出动
        if pot(j, m) <= 0             % 物理不可行组合（精确扫描无窗口）
            ok = false; break;
        end
        covered(m) = true;
        sc = sc + pot(j, m);
    end
    if ok && all(covered)             % 每枚导弹至少被 1 机干扰
        valid(c) = true;
        scores(c) = sc;
        assigns{c} = a;
    end
end

% ---- 按目标值排序，取最优 + 次优 ----
idx = find(valid);
[~, ord] = sort(scores(idx), 'descend');
ord = idx(ord);
nGet = min(k, numel(ord));

for t = 1:nGet
    a = assigns{ord(t)};
    for m = 1:nM
        srv = zeros(0, 2);
        for j = 1:nU
            if a(j) == m
                srv(end+1, :) = [j, 3];   % 每机允许投满 3 弹（第2层自动筛选有效弹位）
            end
        end
        cand(t).missile{m} = srv;
    end
end

%% ---- 控制台输出候选方案 ----
uav_names = {'FY1','FY2','FY3','FY4','FY5'};
for t = 1:nGet
    fprintf('--- 候选指派 %d（ILP 目标 Σpot=%.2f）---\n', t, scores(ord(t)));
    for m = 1:nM
        srv = cand(t).missile{m};
        if isempty(srv)
            fprintf('  M%d <- (无服务机)\n', m);
        else
            fprintf('  M%d <- ', m);
            for r = 1:size(srv, 1)
                fprintf('%s x %d弹  ', uav_names{srv(r, 1)}, srv(r, 2));
            end
            fprintf('\n');
        end
    end
end
end
