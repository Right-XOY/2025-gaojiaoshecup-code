%% 3D投影法有效遮挡示意图｜左→右：导弹‑烟幕球‑圆柱目标
%% 圆柱【左半边侧面白色(可被照到)，右半边侧面黑色(照不到)】
clear; clc; close all;

%% 几何参数：从左到右排布 P1(导弹) --- C(烟幕球) --- P2(圆柱目标)
P1 = [-10, 0, 0];        % 导弹P1【最左侧】
C  = [0, 0, 0];          % 烟幕球心C【中间】
r_smoke = 2.5;           % 烟幕球体半径
P2 = [8, -1.5, 0];       % 真目标圆柱中心P2【最右侧】
cyl_r = 1.4;             % 圆柱半径
cyl_h = 2.2;             % 圆柱高度

fig = figure('Color','w');
ax = axes('Parent',fig);
hold(ax,'on');
axis equal;

%% 光源
light('Position',[-12 8 10],'Style','infinite');
lighting gouraud;
material dull;

%% 1.绘制烟幕球体
[Xsp,Ysp,Zsp] = sphere(60);
Xsp = Xsp*r_smoke + C(1);
Ysp = Ysp*r_smoke + C(2);
Zsp = Zsp*r_smoke + C(3);
surf(ax,Xsp,Ysp,Zsp,'FaceColor',[0.75,0.75,0.75],'FaceAlpha',0.7,'EdgeColor','none');

% 球心C、半径r标注
scatter3(ax,C(1),C(2),C(3),40,'k','filled');
text(ax,C(1)+0.4,C(2)+0.3,C(3)+0.3,'\itC','FontSize',13,'FontName','Times New Roman');
plot3(ax,[C(1),C(1)+r_smoke],[C(2),C(2)],[C(3),C(3)],'k-','LineWidth',1.2);
text(ax,C(1)+1.0,C(2)+0.2,C(3)+0.3,'\itr','FontSize',13,'FontName','Times New Roman');

%% 2.导弹P1（最左边）
scatter3(ax,P1(1),P1(2),P1(3),50,'k','filled');
text(ax,P1(1)-0.9,P1(2)+0.4,P1(3)+0.2,'\itP_1','FontSize',14,'FontName','Times New Roman');

%% 3.计算两条外切线（P1到烟幕球的切线）
dx = C(1)-P1(1);
dy = C(2)-P1(2);
dz = C(3)-P1(3);
d = sqrt(dx^2+dy^2+dz^2);
alpha = asin(r_smoke/d);
beta = atan2(dy,dx);

L = 18;
ang1 = beta + alpha;
ang2 = beta - alpha;
t1x = P1(1)+L*cos(ang1); t1y = P1(2)+L*sin(ang1); t1z = P1(3);
t2x = P1(1)+L*cos(ang2); t2y = P1(2)+L*sin(ang2); t2z = P1(3);

plot3(ax,[P1(1),t1x],[P1(2),t1y],[P1(3),t1z],'k-','LineWidth',1.2);
plot3(ax,[P1(1),t2x],[P1(2),t2y],[P1(3),t2z],'k-','LineWidth',1.2);

%% 4.主视线 P1 -> P2 虚线
plot3(ax,[P1(1),P2(1)],[P1(2),P2(2)],[P1(3),P2(3)],'k--','LineWidth',1.2);

%% 5.绘制真目标圆柱
theta = linspace(0,2*pi,60);
x_cyl = cyl_r*cos(theta) + P2(1);
y_cyl = cyl_r*sin(theta) + P2(2);
z_low  = P2(3) - cyl_h/2;
z_high = P2(3) + cyl_h/2;

% =========圆柱左半边侧面：白色，可被照到=========
theta_illum = linspace(pi/2, 3*pi/2, 80);
xs_illum = P2(1) + cyl_r*cos(theta_illum);
ys_illum = P2(2) + cyl_r*sin(theta_illum);
X_illum = [xs_illum; xs_illum];
Y_illum = [ys_illum; ys_illum];
Z_illum = [z_low * ones(size(xs_illum)); z_high * ones(size(xs_illum))];
surf(ax, X_illum, Y_illum, Z_illum, 'FaceColor', 'w', 'EdgeColor', 'k');

% =========圆柱右半边侧面：黑色，照不到=========
theta_shadow = linspace(-pi/2, pi/2, 80);
xs_shadow = P2(1) + cyl_r*cos(theta_shadow);
ys_shadow = P2(2) + cyl_r*sin(theta_shadow);
X_shadow = [xs_shadow; xs_shadow];
Y_shadow = [ys_shadow; ys_shadow];
Z_shadow = [z_low * ones(size(xs_shadow)); z_high * ones(size(xs_shadow))];
surf(ax, X_shadow, Y_shadow, Z_shadow, 'FaceColor', 'k', 'EdgeColor', 'k');

% 圆柱上下底面
fill3(ax, x_cyl, y_cyl, z_low * ones(size(x_cyl)), 'w', 'EdgeColor', 'k');
fill3(ax, x_cyl, y_cyl, z_high * ones(size(x_cyl)), 'w', 'EdgeColor', 'k');

scatter3(ax, P2(1), P2(2), P2(3), 40, 'k', 'filled');
text(ax, P2(1)+0.5, P2(2)-0.6, P2(3), '\itP_2','FontSize',14,'FontName','Times New Roman');

% 外侧标注
text(ax, P2(1)+2.3, P2(2)-2.3, P2(3), '黑色区域：照不到','FontSize',11,'FontName','SimHei');

%% 视角：稍微俯视，物体从左到右排布
view(ax,-25,30);
grid off;
box off;
xticks([]); yticks([]); zticks([]);

%% 输出，传入fig句柄修复print报错
print(fig,'proj_3d_occlusion','-depsc');
savefig(fig,'proj_3d_occlusion.fig');

hold(ax,'off');
