%% q3_3d.m  问题3 三维策略态势图（风格同 background_scene_overview.png）
%  FY1 投放 3 枚干扰弹干扰 M1。
clear; clc; close all;
cd(fileparts(mfilename('fullpath')));

%% ================== 1. 场景常量 ==================
M0 = [20000; 0; 2000];       % M1
U0 = [17800; 0; 1800];       % FY1
fake   = [0; 0; 0];
true_c = [0; 200; 0];
Rc = 7;  Hc = 10;
v_m = 300;

c_m    = [0.85 0.20 0.20];
c_u    = [0.15 0.40 0.75];
c_fake = [0.15 0.40 0.75];
c_true = [0.15 0.40 0.75];
c_drop = [0.30 0.65 0.30];
c_det  = [0.15 0.15 0.15];
c_cloud= [0.60 0.60 0.60];
fs = 10;

%% ================== 2. 问题3 最优策略（FY1 三弹） ==================
theta = 179.6393;  v = 137.9157;
u_h = [cosd(theta); sind(theta); 0];

t_l = [0.0545, 3.2598, 5.3211];
tau = [3.6396, 5.1123, 5.9222];
t_d = [3.6940, 8.3721, 11.2433];
P_l = [17792.49  17350.44  17066.15;
          0.05      2.83      4.62;
       1800.00   1800.00   1800.00];
P_d = [17290.55  16645.38  16249.40;
          3.21      7.27      9.76;
       1735.09   1671.93   1628.14];

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

% 三枚弹：示意烟幕球（实际半径 10 m，此处示意放大为 250 m）、云团下沉
% daspect([6 2 1]) 会将 X 拉伸 6 倍、Y 拉伸 2 倍，
% 数据坐标乘以 [6 2 1] 可补偿为视觉正球
for i = 1:3
    [Xs, Ys, Zs] = sphere(30);
    surf(P_d(1,i) + 250*6*Xs, P_d(2,i) + 250*2*Ys, P_d(3,i) + 250*Zs, ...
         'FaceColor', c_cloud, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
    cdt = linspace(t_d(i), min(t_d(i)+20, T_end), 80);
    Ctraj = P_d(:,i) - [0;0;3]*(cdt - t_d(i));
    plot3(Ctraj(1,:), Ctraj(2,:), Ctraj(3,:), '--', 'Color', c_cloud, 'LineWidth', 1.2);
end

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


% 坐标轴刻度
set(gca, 'FontName', 'Times New Roman', 'FontSize', fs, 'LineWidth', 0.8, ...
    'XLim', [0 20000], 'YLim', [-3000 2000], 'ZLim', [0 2000], ...
    'XTick', 0:5000:20000, 'XTickLabel', {'0','5000','10000','15000','20000'}, ...
    'YTick', -3000:1000:2000, 'YTickLabel', {'-3000','-2000','-1000','0','1000','2000'}, ...
    'ZTick', 0:500:2000, 'ZTickLabel', {'0','500','1000','1500','2000'});
box on;

print(gcf, 'q3_3d', '-dpng', '-r300');
fprintf('已输出：q3_3d.png\n');
