%% sensitivity_n_convergence.m
%% X对数轴：自动放大n=6,12,24,48；后端100-600压缩；无网格
clear; clc; close all;

n_list  = [6, 12, 24, 48, 100, 200, 300, 400, 600];
T_rec   = [1.392644, 1.391983, 1.391687, 1.391677, 1.391643, ...
           1.391643, 1.391643, 1.391643, 1.391643];
T_inf   = 1.391643;

figure('Color','w','Position',[100,100,760,420]);

fs = 10.5;
lw_curve = 1.6;
lw_ref   = 1.2;
ms       = 6.5;

semilogx(n_list, T_rec, 'o-', ...
    'Color', [0.10, 0.10, 0.10], ...
    'MarkerFaceColor', [0.10, 0.10, 0.10], ...
    'MarkerEdgeColor', 'w', ...
    'MarkerSize', ms, ...
    'LineWidth', lw_curve);
hold on;

% 收敛基准线
yline(T_inf, '--', ...
    'Color', [0.65, 0.25, 0.25], ...
    'LineWidth', lw_ref);

%% 坐标轴设置
xlabel('圆周采样点数 \it n', 'FontSize', fs);
ylabel('有效遮蔽时长 \it T \rm (s)', 'FontSize', fs);

xlim([5, 650]);
ylim([1.3910, 1.3930]);

% 强制显示全部原始n作为刻度
xticks(n_list);
xticklabels(string(n_list));

grid off;       % 彻底关闭网格
set(gca, 'LineWidth',0.85);
set(gca, 'FontName','Times New Roman');
set(gca, 'FontSize',fs);
set(gca, 'Box','on');

legend({'计算遮蔽时长','收敛极限 \it T\rm_\infty = 1.391643 s'},...
    'Location','southwest','Box','off','FontSize',9);

print(gcf,'sensitivity_n_convergence','-dpng','-r300');
fprintf('输出完成，对数X轴，前几个点自动拉开间距\n');
