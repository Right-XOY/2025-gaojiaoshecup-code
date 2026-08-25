function [x_best, f_best, hist_best, hist_muF, hist_muCR] = q3_jade(fun, lb, ub, opts)
% q3_jade  JADE：带外部存档与自适应控制参数的差分进化算法
%   参考文献：Zhang, J., & Sanderson, A. C. (2009). JADE: Adaptive
%   differential evolution with optional external archive. IEEE TEVC.
%
%   核心机制：
%     1) DE/current-to-pbest/1 变异，引导向量取自当代前 p% 优秀个体；
%     2) 外部存档存储被淘汰个体，变异差分项 x_r1 - x_r2 中的 x_r2
%        可从"种群 ∪ 存档"中选取，增强多样性、避免过早收敛；
%     3) 缩放因子 F 由柯西分布生成，交叉率 CR 由正态分布生成，
%        每代按成功个体的均值（F 用 Lehmer 均值）自适应更新 μ_F、μ_CR；
%     4) 二项交叉 + 贪婪选择（带精英保留效果）。
%
% 输入：
%   fun   - 目标函数句柄，fun(x) 返回标量（最小化）
%   lb,ub - 决策变量下/上界（1 x D）
%   opts  - 结构体：NP（种群规模）、MAXITER（最大代数）、
%           p（p-best 比例，默认 0.05）、c（学习速率，默认 0.1）、
%           X0（可选初始种群，前若干行会被采用）
% 输出：
%   x_best   - 最优解
%   f_best   - 最优目标值
%   hist_best - 每代最优目标值（MAXITER x 1）
%   hist_muF / hist_muCR - 每代自适应均值

NP      = opts.NP;
MAXITER = opts.MAXITER;
if isfield(opts, 'p'), p = opts.p; else, p = 0.05; end
if isfield(opts, 'c'), c = opts.c; else, c = 0.1;  end
D = numel(lb);

%% ---- 初始化 ----
pop = lb + (ub - lb).*rand(NP, D);
if isfield(opts, 'X0') && ~isempty(opts.X0)
    n0 = min(size(opts.X0, 1), NP);
    pop(1:n0, :) = opts.X0(1:n0, :);
end
fpop = zeros(NP, 1);
for i = 1:NP
    fpop(i) = fun(pop(i, :));
end

archive = zeros(0, D);          % 外部存档（被淘汰个体）
mu_F  = 0.5;                    % 缩放因子均值
mu_CR = 0.5;                    % 交叉率均值

hist_best  = nan(MAXITER, 1);
hist_muF   = nan(MAXITER, 1);
hist_muCR  = nan(MAXITER, 1);

for g = 1:MAXITER
    [~, order] = sort(fpop);
    pnum = max(2, round(p*NP));                 % p-best 个体数
    S_F  = [];  S_CR = [];                      % 成功参数集合

    for i = 1:NP
        % (1) 自适应控制参数
        CR_i = mu_CR + 0.1*randn;
        CR_i = max(0, min(1, CR_i));
        F_i = mu_F + 0.1*tan(pi*(rand - 0.5));  % 柯西分布
        while F_i <= 0
            F_i = mu_F + 0.1*tan(pi*(rand - 0.5));
        end
        F_i = min(F_i, 1);

        % (2) 选择 p-best 引导向量
        x_i  = pop(i, :);
        x_pb = pop(order(randi(pnum)), :);

        % (3) 选择两个互异且不同于 i 的随机个体（r2 可来自存档）
        r1 = randi(NP);
        while r1 == i, r1 = randi(NP); end
        while true
            r2 = randi(NP + size(archive, 1));
            if r2 <= NP && (r2 == i || r2 == r1), continue; end
            break;
        end
        if r2 <= NP
            xr2 = pop(r2, :);
        else
            xr2 = archive(r2 - NP, :);
        end

        % (4) 变异：DE/current-to-pbest/1（含存档项）
        v = x_i + F_i*(x_pb - x_i) + F_i*(pop(r1, :) - xr2);

        % (5) 边界处理：越界分量随机重置到可行域内
        out = (v < lb) | (v > ub);
        if any(out)
            n_out = sum(out);
            v(out) = lb(out) + (ub(out) - lb(out)).*rand(1, n_out);
        end

        % (6) 二项交叉
        j_rand = randi(D);
        mask = rand(1, D) < CR_i;
        mask(j_rand) = true;
        u = x_i;
        u(mask) = v(mask);

        % (7) 贪婪选择 + 存档更新
        fu = fun(u);
        if fu <= fpop(i)
            archive = [archive; x_i];           % 被淘汰个体入存档
            if size(archive, 1) > NP            % 存档上限 NP，随机删减
                rdel = randperm(size(archive, 1), size(archive, 1) - NP);
                archive(rdel, :) = [];
            end
            pop(i, :) = u;
            fpop(i)   = fu;
            S_F  = [S_F, F_i];
            S_CR = [S_CR, CR_i];
        end
    end

    % (8) 自适应更新均值
    if ~isempty(S_CR)
        mu_CR = (1 - c)*mu_CR + c*mean(S_CR);
    end
    if ~isempty(S_F)
        mu_F = (1 - c)*mu_F + c*sum(S_F.^2)/sum(S_F);   % Lehmer 均值
    end

    [f_best, bi] = min(fpop);
    hist_best(g)  = f_best;
    hist_muF(g)   = mu_F;
    hist_muCR(g)  = mu_CR;
end

[f_best, bi] = min(fpop);
x_best = pop(bi, :);
end
