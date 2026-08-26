%% q4_3d.m  问题4 三维策略态势图（风格同 background_scene_overview.png）
%  配色：导弹=红，无人机/假目标/真目标=蓝；轨迹线均用虚线；
%  缩放比例 daspect([6 2 1])，视角 view([-38 26])，Times New Roman。
%  内容：FY1/FY2/FY3 各投放 1 弹接力遮蔽 M1。
clear; clc; close all;
cd(fileparts(mfilename('fullpath')));

%% ================== 1. 场景常量 ==================
M1 = [20000; 0; 2000];       % M1 初始位置
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

%% ================== 2. 问题4 最优策略（三机接力） ==================
u_name = {'FY1','FY2','FY3'};
% 各无人机航向(deg)、速度、投放时刻 t_l、起爆时刻 t_d
u_theta = [  4.60; 226.18;  90.00];
u_v     = [ 79.82; 117.64; 110.00];
u_tl    = [ 0.993;  7.970;  23.968];
u_td    = [ 1.439; 15.977;  27.968];

% 无人机初始位置（每列一架：FY1/FY2/FY3）
U0 = [17800 12000  6000;
         0  1400 -3000;
      1800  1400   700];

% 投放点 P_l、起爆点 P_d（每列一架）
P_l = [17879.01  11350.82  6000.00;
          6.36    723.51  -363.52;
       1800.00   1400.00   700.00];
P_d = [17914.49  10698.62  6000.00;
          9.21     43.88    76.48;
       1799.03   1085.85   621.60];

u_m   = -M1 / norm(M1);
T_end = norm(M1) / v_m;

%% ================== 3. 绘图 ==================
figure('Color','w','Position',[60 60 720 500]);
hold on; grid on; view([-38, 26]);
daspect([6 2 1]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

% 导弹轨迹（红虚线）
mt = linspace(0, T_end, 300);
Mtraj = M1' + v_m * mt' * u_m';
plot3(Mtraj(:,1), Mtraj(:,2), Mtraj(:,3), '--', 'Color', c_m, 'LineWidth', 1.2);
scatter3(M1(1), M1(2), M1(3), 95, c_m, '^', 'filled');
text(M1(1), M1(2), M1(3)+200, 'M1', 'Color', c_m, 'FontSize', fs, 'FontWeight','bold');

% 各无人机：轨迹（蓝虚线，固定角度/速度持续前飞）+ 起爆点 + 云团下沉
for i = 1:3
    u_h = [cosd(u_theta(i)); sind(u_theta(i)); 0];

    % 无人机轨迹（从初始位置持续前飞）
    ut = linspace(0, 200, 100);
    UAVtraj = U0(:,i)' + u_v(i) * ut' * u_h';
    plot3(UAVtraj(:,1), UAVtraj(:,2), UAVtraj(:,3), '--', 'Color', c_u, 'LineWidth', 1.2);
    scatter3(U0(1,i), U0(2,i), U0(3,i), 70, c_u, 'o', 'filled');
    text(U0(1,i), U0(2,i), U0(3,i)+200, u_name{i}, 'Color', c_u, 'FontSize', fs);

    % 示意烟幕球（实际半径 10 m，此处示意放大为 250 m）
    % daspect([6 2 1]) 会将 X 拉伸 6 倍、Y 拉伸 2 倍，
    % 数据坐标乘以 [6 2 1] 可补偿为视觉正球
    [Xs, Ys, Zs] = sphere(30);
    surf(P_d(1,i) + 250*6*Xs, P_d(2,i) + 250*2*Ys, P_d(3,i) + 250*Zs, ...
         'FaceColor', c_cloud, 'FaceAlpha', 0.35, 'EdgeColor', 'none');

    % 云团下沉轨迹（灰虚线，起爆后 3 m/s 下沉）
    cdt = linspace(u_td(i), min(u_td(i)+20, T_end), 80);
    Ctraj = P_d(:,i) - [0;0;3]*(cdt - u_td(i));
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


% 坐标轴刻度（消去科学计数）
set(gca, 'FontName', 'Times New Roman', 'FontSize', fs, 'LineWidth', 0.8, ...
    'XLim', [0 20000], 'YLim', [-3000 2000], 'ZLim', [0 2000], ...
    'XTick', 0:5000:20000, 'XTickLabel', {'0','5000','10000','15000','20000'}, ...
    'YTick', -3000:1000:2000, 'YTickLabel', {'-3000','-2000','-1000','0','1000','2000'}, ...
    'ZTick', 0:500:2000, 'ZTickLabel', {'0','500','1000','1500','2000'});
box on;

print(gcf, 'q4_3d', '-dpng', '-r300');
fprintf('已输出：q4_3d.png\n');
