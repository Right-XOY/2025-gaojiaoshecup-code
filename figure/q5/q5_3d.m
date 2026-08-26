%% q5_3d.m  问题5 三维策略态势图（风格同 background_scene_overview.png）
%  配色：导弹=红，无人机/假目标/真目标=蓝；轨迹线均用虚线；
%  缩放比例 daspect([6 2 1])，视角 view([-38 26])，Times New Roman。
%  内容：5 架无人机（每机至多 3 弹）干扰 M1/M2/M3。
clear; clc; close all;
cd(fileparts(mfilename('fullpath')));

%% ================== 1. 场景常量 ==================
% 来袭导弹初始位置（每列一枚：M1/M2/M3）
M0 = [20000 19000 18000;
         0   600  -600;
      2000  2100  1900];
m_name = {'M1','M2','M3'};

% 无人机初始位置（每列一架：FY1~FY5）
U0 = [17800 12000  6000 11000 13000;
         0  1400 -3000  2000 -2000;
      1800  1400   700  1800  1300];
u_name = {'FY1','FY2','FY3','FY4','FY5'};

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

%% ================== 2. 问题5 最优策略 ==================
% 各无人机航向(deg)、速度、飞行时长（至最后一弹投放）
u_theta = [ 10.0; 270.0;  90.0; 280.0; 120.0];
u_v     = [ 98.85; 128.5; 120.0; 130.0; 115.69];
u_tmax  = [  1.0;   6.0;  22.0;   2.0;  15.0];

% 各无人机投放点（有效弹位，绿菱形），cell 每列一弹 [x;y;z]
drop = cell(1,5);
drop{1} = [17800.0  17897.3;  0.0  17.2;  1800.0  1800.0];           % FY1
drop{2} = [12000.0  12000.0  12000.0; 886.0  757.5  629.0; 1400.0  1400.0  1400.0]; % FY2
drop{3} = [ 6000.0; -360.0;  700.0];                                 % FY3
drop{4} = [11045.1; 1743.9; 1800.0];                                 % FY4
drop{5} = [12190.2  12132.3; -597.4  -497.2; 1300.0  1300.0];        % FY5

% 起爆点/起爆时刻（有效弹），按导弹分组；每列一枚弹：[t_d; x; y; z]
det_M1 = [  0.144     1.044    25.5;
        17814.0   17901.6  6000.0;
            2.5      17.9    60.0;
         1799.9    1800.0   640.0];
det_M2 = [  7.566     7.813     7.963    12.5;
        12000.0   12000.0   12000.0 11282.2;
          427.7     395.9     376.7   399.7;
         1337.7    1361.2    1381.1  1259.8];
det_M3 = [ 16.042    15.931;
        12072.1   12078.5;
         -392.8    -403.9;
         1279.6    1295.7];
det_all = {det_M1, det_M2, det_M3};

% 各导弹到达假目标时刻
T_end = zeros(1,3);
for i = 1:3
    T_end(i) = norm(M0(:,i)) / v_m;
end

%% ================== 3. 绘图 ==================
figure('Color','w','Position',[60 60 720 500]);
hold on; grid on; view([-38, 26]);
daspect([6 2 1]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

% 三枚导弹轨迹（红虚线）
for i = 1:3
    u_m = -M0(:,i) / norm(M0(:,i));
    mt = linspace(0, T_end(i), 300);
    Mtraj = M0(:,i)' + v_m * mt' * u_m';
    plot3(Mtraj(:,1), Mtraj(:,2), Mtraj(:,3), '--', 'Color', c_m, 'LineWidth', 1.2);
    scatter3(M0(1,i), M0(2,i), M0(3,i), 95, c_m, '^', 'filled');
    text(M0(1,i), M0(2,i), M0(3,i)+200, m_name{i}, 'Color', c_m, ...
         'FontSize', fs, 'FontWeight','bold');
end

% 各无人机轨迹（蓝虚线，固定角度/速度持续前飞）
for i = 1:5
    u_h = [cosd(u_theta(i)); sind(u_theta(i)); 0];
    ut = linspace(0, 200, 100);
    UAVtraj = U0(:,i)' + u_v(i) * ut' * u_h';
    plot3(UAVtraj(:,1), UAVtraj(:,2), UAVtraj(:,3), '--', 'Color', c_u, 'LineWidth', 1.2);
    scatter3(U0(1,i), U0(2,i), U0(3,i), 70, c_u, 'o', 'filled');
    text(U0(1,i), U0(2,i), U0(3,i)+200, u_name{i}, 'Color', c_u, 'FontSize', fs);
end

% 各起爆点：示意烟幕球（实际半径 10 m，此处示意放大为 250 m）+ 云团下沉（灰虚线），按导弹分组
% daspect([6 2 1]) 会将 X 拉伸 6 倍、Y 拉伸 2 倍，
% 数据坐标乘以 [6 2 1] 可补偿为视觉正球
for i = 1:3
    det = det_all{i};
    for j = 1:size(det, 2)
        t_d = det(1,j);
        P_d = det(2:4,j);
        [Xs, Ys, Zs] = sphere(30);
        surf(P_d(1) + 250*6*Xs, P_d(2) + 250*2*Ys, P_d(3) + 250*Zs, ...
             'FaceColor', c_cloud, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        cdt = linspace(t_d, min(t_d+20, T_end(i)), 80);
        Ctraj = P_d - [0;0;3]*(cdt - t_d);
        plot3(Ctraj(1,:), Ctraj(2,:), Ctraj(3,:), '--', 'Color', c_cloud, 'LineWidth', 1.2);
    end
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

print(gcf, 'q5_3d', '-dpng', '-r300');
fprintf('已输出：q5_3d.png\n');
