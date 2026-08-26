# -*- coding: utf-8 -*-
"""occlusion_core.py
2025 高教社杯 A 题「烟幕干扰弹的投放策略」论文配图共享模块。

封装：
  1) 场景常量（导弹/无人机初始位置、真目标圆柱、烟幕参数）；
  2) 运动学模型（导弹位置、无人机投放/起爆点、云团中心）；
  3) 遮蔽判定模型（点遮蔽 -> 圆周全遮蔽 -> 连续裕量 H(t)）；
  4) 时间线并集计算（多枚弹对同一导弹的遮蔽并集）。

遮蔽判定准则（与 code/q2_occlusion.m、q3_occlusion_window.m 完全一致）：
    圆柱被完全遮蔽  <=>  上、下底面圆周的所有采样点均被烟幕球遮挡；
    点 P 被遮挡      <=>  烟幕球（中心 C、半径 R）与线段 MP 相交，
                          即 C 到线段 MP 的最近距离 <= R（垂足参数截断到 [0,1]）。
"""

import numpy as np

# ---------------- 场景常量 ----------------
G = 9.8                 # 重力加速度 (m/s^2)
R_SMOKE = 10.0          # 烟幕有效遮蔽半径 (m)
V_SINK = 3.0            # 云团下沉速度 (m/s)
T_SMOKE = 20.0          # 起爆后有效时长 (s)
V_MISSILE = 300.0       # 导弹速度 (m/s)

CYL_CENTER = np.array([0.0, 200.0])   # 真目标圆柱轴心 (x, y)
CYL_R = 7.0                            # 圆柱半径 (m)
CYL_TOP = 10.0                         # 圆柱顶面高度 (m)
CYL_BOTTOM = 0.0                       # 圆柱底面高度 (m)
N_PHI = 100                            # 圆周采样点数（n>=100 已收敛）

# 导弹初始位置 (m)
M0 = {
    'M1': np.array([20000.0, 0.0, 2000.0]),
    'M2': np.array([19000.0, 600.0, 2100.0]),
    'M3': np.array([18000.0, -600.0, 1900.0]),
}

# 无人机初始位置 (m)
UAV0 = {
    'FY1': np.array([17800.0, 0.0, 1800.0]),
    'FY2': np.array([12000.0, 1400.0, 1400.0]),
    'FY3': np.array([6000.0, -3000.0, 700.0]),
    'FY4': np.array([11000.0, 2000.0, 1800.0]),
    'FY5': np.array([13000.0, -2000.0, 1300.0]),
}


# ---------------- 运动学 ----------------
def missile_dir(m0):
    """导弹飞行方向单位向量（直指假目标原点）。"""
    return -m0 / np.linalg.norm(m0)


def missile_pos(m0, t):
    """导弹在 t 时刻的位置。"""
    return m0 + V_MISSILE * t * missile_dir(m0)


def t_end(m0):
    """导弹到达假目标（原点）的时刻。"""
    return np.linalg.norm(m0) / V_MISSILE


def uav_drop_point(uav0, theta_deg, v, t_l):
    """无人机投放点：等高度匀速直线飞行 t_l 秒后的位置。"""
    th = np.deg2rad(theta_deg)
    u_h = np.array([np.cos(th), np.sin(th), 0.0])
    return uav0 + v * t_l * u_h


def uav_detonate_point(uav0, theta_deg, v, t_l, tau):
    """干扰弹起爆点：投放后仅受重力（竖直初速为 0），tau 秒后起爆。"""
    th = np.deg2rad(theta_deg)
    u_h = np.array([np.cos(th), np.sin(th), 0.0])
    P_l = uav0 + v * t_l * u_h
    return P_l + v * tau * u_h - np.array([0.0, 0.0, 0.5 * G * tau ** 2])


def cloud_center(P_d, t, t_d):
    """云团中心：起爆后以 3 m/s 匀速下沉。"""
    return P_d - np.array([0.0, 0.0, V_SINK * (t - t_d)])


# ---------------- 遮蔽判定 ----------------
def _rim_points():
    """真目标上、下底面圆周采样点，返回 3 x (2*N_PHI)。"""
    phi = np.linspace(0.0, 2.0 * np.pi, N_PHI, endpoint=False)
    top = np.vstack([CYL_CENTER[0] + CYL_R * np.cos(phi),
                     CYL_CENTER[1] + CYL_R * np.sin(phi),
                     np.full(N_PHI, CYL_TOP)])
    bot = np.vstack([CYL_CENTER[0] + CYL_R * np.cos(phi),
                     CYL_CENTER[1] + CYL_R * np.sin(phi),
                     np.full(N_PHI, CYL_BOTTOM)])
    return np.hstack([top, bot])          # 3 x 2N


def margin(M, C, R=R_SMOKE):
    """连续裕量 H(t) = max_p ( d_p - R )，H<=0 等价于圆柱被完全遮蔽。"""
    P = _rim_points()
    v = P - M[:, None]                    # 3 x N，导弹 -> 圆周点
    w = C - M                             # 3，导弹 -> 云团中心
    dvv = np.sum(v * v, axis=0)
    tp = (w @ v) / dvv
    tp = np.clip(tp, 0.0, 1.0)            # 垂足参数截断到 [0,1]（到线段距离）
    foot = M[:, None] + tp[None, :] * v
    d = np.sqrt(np.sum((C[:, None] - foot) ** 2, axis=0))
    return float(np.max(d - R))


def occluded(M, C, R=R_SMOKE):
    """布尔遮蔽判定：圆柱是否被完全遮蔽。"""
    return margin(M, C, R) <= 0.0


# ---------------- 时间线 ----------------
def occlusion_timeline(m0, shells, dt=0.005):
    """多枚弹对同一导弹的遮蔽并集时间线。

    shells : list of (t_d, P_d)，每枚弹的起爆时刻与起爆点。
    返回 (ts, shielded)，ts 为 0..t_end 时间网格，shielded 为布尔向量。
    """
    te = t_end(m0)
    ts = np.arange(0.0, te + dt, dt)
    shielded = np.zeros(ts.shape, dtype=bool)
    for (t_d, P_d) in shells:
        lo = int(np.floor(t_d / dt))
        hi = int(np.ceil((t_d + T_SMOKE) / dt))
        lo = max(0, lo)
        hi = min(len(ts), hi + 1)
        for k in range(lo, hi):
            t = ts[k]
            if t < t_d or t > t_d + T_SMOKE:
                continue
            if occluded(missile_pos(m0, t), cloud_center(P_d, t, t_d)):
                shielded[k] = True
    return ts, shielded


def occlusion_intervals(ts, shielded, dt):
    """从布尔时间线提取遮蔽区间 [start, stop]。"""
    edges = np.diff(np.concatenate([[0], shielded.astype(int), [0]]))
    starts = ts[np.where(edges == 1)[0]]
    stops = ts[np.where(edges == -1)[0] - 1]
    return list(zip(starts, stops))
