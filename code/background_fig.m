%% 三维场景绘图：空地导弹、无人机、真假目标
clear; clc; close all;

%% 1. 坐标数据 (单位：m)
% 导弹 M1,M2,M3
missile = [
    20000,   0, 2000;   % M1
    19000, 600, 2100;   % M2
    18000,-600, 1900];  % M3
m_name = {'M1','M2','M3'};

% 无人机 FY1~FY5
uav = [
    17800,    0, 1800;  % FY1
    12000, 1400, 1400;  % FY2
    6000, -3000,  700;  % FY3
    11000, 2000, 1800;  % FY4
    13000,-2000, 1300]; % FY5
u_name = {'FY1','FY2','FY3','FY4','FY5'};

false_target = [0,0,0];       % 假目标原点
true_target_center = [0,200,0];% 真目标下底面圆心
R = 7;    % 真目标圆柱半径
H = 10;   % 真目标圆柱高度

%% 2. 创建画布
figure('Color','w');
hold on; grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('导弹‑无人机‑真假目标三维态势图');

%% 绘制导弹 红色上三角
scatter3(missile(:,1),missile(:,2),missile(:,3),80,'r','^','filled');
for i = 1:size(missile,1)
    text(missile(i,1)+300, missile(i,2), missile(i,3), m_name{i},'Color','r','FontSize',9);
end

%% 绘制导弹指向假目标原点的飞行射线
for i = 1:size(missile,1)
    pts = [missile(i,:); false_target];
    plot3(pts(:,1),pts(:,2),pts(:,3),'r--','LineWidth',1);
end

%% 绘制无人机 蓝色圆点
scatter3(uav(:,1),uav(:,2),uav(:,3),60,'b','o','filled');
for i = 1:size(uav,1)
    text(uav(i,1)+300, uav(i,2), uav(i,3), u_name{i},'Color','b','FontSize',9);
end

%% 绘制假目标原点 橙色星号
scatter3(false_target(1),false_target(2),false_target(3),100,'m','*');
text(false_target(1)+500,false_target(2),false_target(3),'假目标(0,0,0)','Color','m','FontSize',10);

%% 绘制真目标圆柱体
theta = linspace(0,2*pi,40);
x_cyl = true_target_center(1) + R*cos(theta);
y_cyl = true_target_center(2) + R*sin(theta);
z0 = true_target_center(3);
z1 = z0 + H;

% 圆柱侧面
Xc = repmat(x_cyl,2,1);
Yc = repmat(y_cyl,2,1);
Zc = [z0*ones(size(x_cyl)); z1*ones(size(x_cyl))];
surf(Xc,Yc,Zc,'FaceColor','g','FaceAlpha',0.3,'EdgeColor','g');
text(true_target_center(1),true_target_center(2),z1+20,'真目标圆柱 r=7m h=10m','Color','g','FontSize',9);

%% 图例、视角
legend('导弹','','无人机','假目标','真目标','Location','best');
view(3);
hold off;