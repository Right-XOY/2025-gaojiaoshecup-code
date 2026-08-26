%% q2_3d.m  问题2 三维策略态势图（风格同 background_scene_overview.png）
%  配色：导弹=红，无人机/假目标/真目标=蓝；轨迹线均用虚线；
%  缩放比例 daspect([6 2 1])，视角 view([-38 26])，Times New Roman。
clear; clc; close all;
cd(fileparts(mfilename('fullpath')));

%% ================== 1. 场景常量 ==================
M0 = [20000; 0; 2000];       % M1 初始位置
U0 = [17800; 0; 1800];       % FY1 初始位置
fake   = [0; 0; 0];          % 假目标（原点）
true_c = [0; 200; 0];        % 真目标下底面圆心
Rc = 7;  Hc = 10;            % 真目标圆柱半径、高
v_m = 300;                   % 导弹速度

% 配色（与 background_scene_overview 一致）
c_m    = [0.85 0.20 0.20];   % 导弹：红
c_u    = [0.15 0.40 0.75];   % 无人机：蓝
c_fake = [0.15 0.40 0.75];   % 假目标：蓝
c_true = [0.15 0.40 0.75];   % 真目标：蓝
c_drop = [0.30 0.65 0.30];   % 投放点：绿
c_det  = [0.15 0.15 0.15];   % 起爆点：黑
c_cloud= [0.60 0.60 0.60];   % 云团下沉：灰
fs = 10;                     % 字号

%% ================== 2. 问题2 最优策略 ==================
theta = 8.15;  v = 86.4;  t_l = 0.001;  tau = 0.976;  t_d = 0.977;
u_h = [cosd(theta); sind(theta); 0];
P_l = [17800.09; 0.01; 1800.00];     % 投放点
P_d = [17883.56; 11.96; 1795.33];    % 起爆点

% 导弹轨迹（直指假目标）
u_m   = -M0 / norm(M0);
T_end = norm(M0) / v_m;

%% ================== 3. 绘图 ==================
figure('Color','w','Position',[60 60 720 500]);
hold on; grid on; view([-38, 26]);
daspect([6 2 1]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

% 导弹轨迹（红虚线）
mt = linspace(0, T_end, 300);
Mtraj = M0' + v_m * mt' * u_m';
plot3(Mtraj(:,1), Mtraj(:,2), Mtraj(:,3), '--', 'Color', c_m, 'LineWidth', 1.2);
scatter3(M0(1), M0(2), M0(3), 95, c_m, '^', 'filled');
text(M0(1), M0(2), M0(3)+200, 'M1', 'Color', c_m, 'FontSize', fs, 'FontWeight','bold');

% 无人机轨迹（蓝虚线，固定角度/速度持续前飞）
ut = linspace(0, 200, 100);
UAVtraj = U0' + v * ut' * u_h';
plot3(UAVtraj(:,1), UAVtraj(:,2), UAVtraj(:,3), '--', 'Color', c_u, 'LineWidth', 1.2);
scatter3(U0(1), U0(2), U0(3), 70, c_u, 'o', 'filled');
text(U0(1), U0(2), U0(3)+200, 'FY1', 'Color', c_u, 'FontSize', fs);

% 起爆点：示意烟幕球（实际半径 10 m，此处示意放大为 250 m）
% daspect([6 2 1]) 会将 X 拉伸 6 倍、Y 拉伸 2 倍，
% 数据坐标乘以 [6 2 1] 可补偿为视觉正球
[Xs, Ys, Zs] = sphere(30);
surf(P_d(1) + 250*6*Xs, P_d(2) + 250*2*Ys, P_d(3) + 250*Zs, ...
     'FaceColor', c_cloud, 'FaceAlpha', 0.35, 'EdgeColor', 'none');

% 云团下沉轨迹（灰虚线）
cdt = linspace(t_d, min(t_d+20, T_end), 80);
Ctraj = P_d - [0;0;3]*(cdt - t_d);
plot3(Ctraj(1,:), Ctraj(2,:), Ctraj(3,:), '--', 'Color', c_cloud, 'LineWidth', 1.2);

% 假目标（蓝圆点）
scatter3(0, 0, 0, 120, c_fake, 'o', 'filled');
text(3000, 0, -220, '假目标', 'Color', c_fake, 'FontSize', fs, 'HorizontalAlignment','center');

% 真目标圆柱（蓝半透明）
theta_c = linspace(0, 2*pi, 80);
x_cyl = true_c(1) + Rc*cos(theta_c);
y_cyl = true_c(2) + Rc*sin(theta_c);
z_bot = true_c(3);  z_top = true_c(3) + Hc;
Xc = repmat(x_cyl, 2, 1);
Yc = repmat(y_cyl, 2, 1);
Zc = [z_bot*ones(size(x_cyl)); z_top*ones(size(x_cyl))];
surf(Xc, Yc, Zc, 'FaceColor', c_true, 'FaceAlpha', 0.45, 'EdgeColor', 'none');
fill3(x_cyl, y_cyl, z_top*ones(size(x_cyl)), c_true, 'FaceAlpha', 0.45, 'EdgeColor','none');
fill3(x_cyl, y_cyl, z_bot*ones(size(x_cyl)), c_true, 'FaceAlpha', 0.45, 'EdgeColor','none');


% 坐标轴刻度（消去科学计数）
set(gca, 'FontName', 'Times New Roman', 'FontSize', fs, 'LineWidth', 0.8, ...
    'XLim', [0 20000], 'YLim', [-3000 2000], 'ZLim', [0 2000], ...
    'XTick', 0:5000:20000, 'XTickLabel', {'0','5000','10000','15000','20000'}, ...
    'YTick', -3000:1000:2000, 'YTickLabel', {'-3000','-2000','-1000','0','1000','2000'}, ...
    'ZTick', 0:500:2000, 'ZTickLabel', {'0','500','1000','1500','2000'});
box on;

print(gcf, 'q2_3d', '-dpng', '-r300');
fprintf('已输出：q2_3d.png\n');
