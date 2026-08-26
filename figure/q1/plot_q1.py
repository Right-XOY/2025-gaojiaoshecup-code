# -*- coding: utf-8 -*-
"""plot_q1.py  问题1 论文配图（4 张）
    1) q1_timeline.png    遮蔽时间线图
    2) q1_margin.png      裕量函数 H(t) 与遮蔽区间
    3) q1_3d.png          三维态势图
    4) q1_projection.png  导弹视角投影示意图
数据来自 结果.md：FY1 以 120 m/s 朝向假目标，1.5 s 投放、3.6 s 后起爆。
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.occlusion_core import (M0, UAV0, missile_pos, cloud_center, t_end,
                            margin, occlusion_timeline, occlusion_intervals,
                            CYL_CENTER, CYL_R, CYL_TOP, CYL_BOTTOM, V_SINK,
                            T_SMOKE, R_SMOKE)

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

HERE = os.path.dirname(os.path.abspath(__file__))


def save(fig, name):
    fig.savefig(os.path.join(HERE, name), dpi=300, bbox_inches='tight')
    plt.close(fig)


def draw_true_target(ax, alpha=0.25):
    """在 3D 轴上画真目标圆柱（半径 7 m、高 10 m，轴心 (0,200,0)）。"""
    theta = np.linspace(0, 2 * np.pi, 60)
    z = np.array([CYL_BOTTOM, CYL_TOP])
    tg, zg = np.meshgrid(theta, z)
    xg = CYL_CENTER[0] + CYL_R * np.cos(tg)
    yg = CYL_CENTER[1] + CYL_R * np.sin(tg)
    ax.plot_surface(xg, yg, zg, color='#2ca02c', alpha=alpha, linewidth=0)
    # 上下底面圆周轮廓
    xc = CYL_CENTER[0] + CYL_R * np.cos(theta)
    yc = CYL_CENTER[1] + CYL_R * np.sin(theta)
    for zz in (CYL_BOTTOM, CYL_TOP):
        ax.plot(xc, yc, np.full_like(xc, zz), color='#2ca02c', linewidth=0.8)


# ---------------- 场景数据 ----------------
M1 = M0['M1']
uav0 = UAV0['FY1']
v_u, t_l, tau = 120.0, 1.5, 3.6
P_l = np.array([17620.0, 0.0, 1800.0])
P_d = np.array([17188.0, 0.0, 1736.496])
t_d = 5.1
t_start, t_stop = 8.056445, 9.448088
T = 1.391643

# ================ 1. 遮蔽时间线 ================
ts, sh = occlusion_timeline(M1, [(t_d, P_d)])
ivs = occlusion_intervals(ts, sh, ts[1] - ts[0])

fig, ax = plt.subplots(figsize=(8, 3.6))
for (a, b) in ivs:
    ax.axvspan(a, b, color='#e64b4b', alpha=0.25, lw=0)
ax.step(ts, sh.astype(int), where='post', color='k', lw=1.5)
ax.axvline(t_d, color='b', ls='--', lw=1, label='起爆时刻')
ax.axvline(t_start, color='m', ls='--', lw=1, label='遮蔽开始')
ax.axvline(t_stop, color='m', ls='--', lw=1, label='遮蔽结束')
ax.set_yticks([0, 1]); ax.set_ylim(-0.15, 1.15)
ax.set_xlabel('时间 t (s)'); ax.set_ylabel('遮蔽状态 I(t)')
ax.set_title(f'M1 对真目标的遮蔽时间线（有效遮蔽时长 = {T:.4f} s）')
ax.grid(True, ls=':', alpha=0.5)
ax.legend(loc='best', fontsize=8)
save(fig, 'q1_timeline.png')

# ================ 2. 裕量函数 H(t) ================
tg = np.linspace(t_d, t_d + T_SMOKE, 800)
H = np.array([margin(missile_pos(M1, t), cloud_center(P_d, t, t_d)) for t in tg])

fig, ax = plt.subplots(figsize=(8, 3.6))
ax.plot(tg, H, 'k-', lw=1.5)
ax.axhline(0, color='r', ls='--', lw=1, label='H=0')
ax.axvline(t_start, color='m', ls='--', lw=1, label='$t_{start}$')
ax.axvline(t_stop, color='m', ls='--', lw=1, label='$t_{end}$')
ax.scatter([t_start, t_stop], [0, 0], s=50, color='r', zorder=5)
ax.set_xlabel('时间 t (s)'); ax.set_ylabel('裕量函数 H(t)')
ax.set_title(f'裕量函数 H(t)（H(t)≤0 即遮蔽；t*=[{t_start:.4f}, {t_stop:.4f}] s）')
ax.grid(True, ls=':', alpha=0.5)
ax.legend(loc='best', fontsize=8)
save(fig, 'q1_margin.png')

# ================ 3. 三维态势图 ================
fig = plt.figure(figsize=(8.4, 6.6))
ax = fig.add_subplot(111, projection='3d')

# 导弹轨迹
mt = np.linspace(0, t_end(M1), 300)
Mt = np.array([missile_pos(M1, t) for t in mt])
ax.plot(Mt[:, 0], Mt[:, 1], Mt[:, 2], 'r-', lw=1.6, label='导弹轨迹')
ax.scatter(*M1, s=45, c='r', marker='^')
ax.text(*M1 + [300, 0, 0], 'M1', color='r', fontsize=10)

# 无人机轨迹（朝向假目标 -x 方向飞 t_l 秒）
u_h = np.array([-1.0, 0.0, 0.0])
ut = np.linspace(0, t_l, 60)
Ut = uav0[None, :] + v_u * ut[:, None] * u_h[None, :]
ax.plot(Ut[:, 0], Ut[:, 1], Ut[:, 2], 'b-', lw=1.6, label='无人机轨迹')
ax.scatter(*uav0, s=40, c='b', marker='o')
ax.text(*uav0 + [300, 0, 0], 'FY1', color='b', fontsize=10)

# 投放点、起爆点
ax.scatter(*P_l, s=45, c='g', marker='D', label='投放点')
ax.scatter(*P_d, s=55, c='k', marker='p', label='起爆点')
ax.text(*P_l + [300, 0, 0], '投放点', color='g', fontsize=9)
ax.text(*P_d + [300, 0, 0], '起爆点', color='k', fontsize=9)

# 云团下沉轨迹
cdt = np.linspace(t_d, t_d + T_SMOKE, 80)
Ct = np.array([cloud_center(P_d, t, t_d) for t in cdt])
ax.plot(Ct[:, 0], Ct[:, 1], Ct[:, 2], 'm-', lw=1.6, label='云团下沉')

# 假目标与真目标
ax.scatter(0, 0, 0, s=45, c='m', marker='*')
ax.text(500, 0, -150, '假目标(0,0,0)', color='m', fontsize=9)
draw_true_target(ax)

ax.set_xlabel('X (m)'); ax.set_ylabel('Y (m)'); ax.set_zlabel('Z (m)')
ax.set_title('问题1 给定投放策略三维态势图')
ax.set_box_aspect([6, 2, 1])
ax.view_init(elev=22, azim=-58)
ax.legend(loc='best', fontsize=7)
save(fig, 'q1_3d.png')

# ================ 4. 导弹视角投影示意图 ================
t_v = 8.7  # 遮蔽窗口内代表时刻
P_m = missile_pos(M1, t_v)
C_v = cloud_center(P_d, t_v, t_d)
T_v = np.array([0.0, 200.0, 0.0])
ww = (T_v - P_m) / np.linalg.norm(T_v - P_m)
uu = np.cross(ww, np.array([0, 0, 1.0]))
uu /= np.linalg.norm(uu)
vv = np.cross(ww, uu)


def projA(X):
    return float((X - P_m) @ uu) / np.linalg.norm(X - P_m)


def projB(X):
    return float((X - P_m) @ vv) / np.linalg.norm(X - P_m)


phi_c = np.linspace(0, 2 * np.pi, 120)
Pt = np.vstack([CYL_CENTER[0] + CYL_R * np.cos(phi_c),
                CYL_CENTER[1] + CYL_R * np.sin(phi_c),
                np.full(120, CYL_TOP)])
Pb = np.vstack([CYL_CENTER[0] + CYL_R * np.cos(phi_c),
                CYL_CENTER[1] + CYL_R * np.sin(phi_c),
                np.full(120, CYL_BOTTOM)])
A_top = np.array([projA(Pt[:, i]) for i in range(120)])
B_top = np.array([projB(Pt[:, i]) for i in range(120)])
A_bot = np.array([projA(Pb[:, i]) for i in range(120)])
B_bot = np.array([projB(Pb[:, i]) for i in range(120)])

A_c, B_c = projA(C_v), projB(C_v)
alpha_v = np.arcsin(R_SMOKE / np.linalg.norm(P_m - C_v))

r_bot = np.sqrt(A_bot ** 2 + B_bot ** 2)
front_arc = np.argsort(r_bot)[:60]
back_arc = np.setdiff1d(np.arange(120), front_arc)

th_d = np.linspace(0, 2 * np.pi, 100)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.6))

ax1.fill(A_c + alpha_v * np.cos(th_d), B_c + alpha_v * np.sin(th_d),
         color='#d5d5d5', alpha=0.5, ec='k', lw=1.1, label='烟幕角圆盘')
ax1.plot(A_top, B_top, 'b-', lw=1.5, label='顶面圆周')
ax1.plot(A_bot[front_arc], B_bot[front_arc], 'g-', lw=1.5, label='底面前弧(可见)')
ax1.plot(A_bot[back_arc], B_bot[back_arc], 'g--', lw=1, label='底面后弧')
ax1.plot([0, A_c], [0, B_c], 'm--', lw=1.2, label='云心方向')
ax1.scatter(0, 0, s=30, c='r', zorder=5)
ax1.scatter(A_c, B_c, s=30, c='m', marker='p', zorder=5)
ax1.set_aspect('equal')
ax1.set_xlabel('水平角偏移 (rad)'); ax1.set_ylabel('垂直角偏移 (rad)')
ax1.set_title(f'导弹视角：目标轮廓与烟幕角圆盘（t={t_v:.1f} s，遮蔽中）')
ax1.grid(True, ls=':', alpha=0.5)
ax1.legend(loc='best', fontsize=7)

ax2.plot(A_top * 1e4, B_top * 1e4, 'b-', lw=1.5, label='顶面圆周')
ax2.plot(A_bot[front_arc] * 1e4, B_bot[front_arc] * 1e4, 'g-', lw=1.5, label='底面前弧')
ax2.scatter(0, 0, s=30, c='r', zorder=5)
ax2.scatter(A_c * 1e4, B_c * 1e4, s=30, c='m', marker='p', zorder=5, label='云心')
ax2.plot([0, A_c * 1e4], [0, B_c * 1e4], 'm--', lw=1.2)
ax2.set_aspect('equal')
ax2.set_xlabel(r'水平角偏移 ($\times 10^{-4}$ rad)')
ax2.set_ylabel(r'垂直角偏移 ($\times 10^{-4}$ rad)')
ax2.set_title('目标轮廓局部放大')
ax2.grid(True, ls=':', alpha=0.5)
ax2.legend(loc='best', fontsize=7)

save(fig, 'q1_projection.png')
print('Q1 四张图已生成：q1_timeline.png / q1_margin.png / q1_3d.png / q1_projection.png')
