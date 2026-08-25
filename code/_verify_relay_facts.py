# -*- coding: utf-8 -*-
"""快速核对：MATLAB 与 Python 遮蔽判定是否一致，接力是否真实。"""
import numpy as np

G = 9.8
R = 10.0
VM = 300.0
nPhi = 100
ang = np.arange(nPhi) * 2 * np.pi / nPhi
rim = []
for z in [0.0, 10.0]:
    for a in ang:
        rim.append([7 * np.cos(a), 200 + 7 * np.sin(a), z])
rim = np.array(rim).T


def shield(uav0, M0, theta_deg, v, t_ls, taus, dt=0.005):
    T_end = np.linalg.norm(M0) / VM
    um = -M0 / np.linalg.norm(M0)
    uh = np.array([np.cos(np.radians(theta_deg)), np.sin(np.radians(theta_deg)), 0.0])
    ts = np.arange(0.0, T_end, dt)
    N = len(ts)
    Ms = M0[:, None] + VM * um[:, None] * ts[None, :]
    shielded = np.zeros(N, dtype=bool)
    for t_l, tau in zip(t_ls, taus):
        Pl = uav0 + v * t_l * uh
        Pd = Pl + v * tau * uh - np.array([0.0, 0.0, 0.5 * G * tau ** 2])
        td = t_l + tau
        C = Pd[:, None] - np.array([0.0, 0.0, 3.0])[:, None] * (ts[None, :] - td)
        vvec = rim[:, None, :] - Ms[:, :, None]
        wvec = C[:, :, None] - Ms[:, :, None]
        vv_norm = np.sum(vvec * vvec, axis=0)
        tp = np.sum(wvec * vvec, axis=0) / vv_norm
        tp = np.clip(tp, 0.0, 1.0)
        closest = Ms[:, :, None] + tp[None, :, :] * vvec
        d = np.sqrt(np.sum((C[:, :, None] - closest) ** 2, axis=0))
        maxd = np.max(d, axis=1)
        shielded |= (ts >= td) & (ts <= td + 20) & (maxd <= R)
    dur = np.sum(shielded) * dt
    if shielded.any():
        idx = np.where(shielded)[0]
        return dur, ts[idx[0]], ts[idx[-1]]
    return dur, None, None


def report(uav0, M0, th, v, t_ls, taus, label):
    d, s, e = shield(uav0, M0, th, v, t_ls, taus)
    print(f"{label}: {d:.4f}s  窗口=[{s if s is None else round(s,3)}, {e if e is None else round(e,3)}]")


# FY1 -> M1
u1 = np.array([17800.0, 0.0, 1800.0])
M1 = np.array([20000.0, 0.0, 2000.0])
print("== FY1->M1 ==")
report(u1, M1, 10, 110, [0], [0.5], "单弹(t_l=0,τ=0.5,v=110)")
report(u1, M1, 10, 110, [0, 1, 2], [0.5, 0.5, 0.5], "三弹gap1同τ(v=110)")
report(u1, M1, 10, 85, [0, 1, 2], [0.3, 0.3, 0.3], "三弹gap1同τ(v=85,τ=0.3)")
report(u1, M1, 180, 110, [0, 1, 2], [0.5, 0.5, 0.5], "朝-x三弹(v=110)")

# FY2 -> M2
u2 = np.array([12000.0, 1400.0, 1400.0])
M2 = np.array([19000.0, 600.0, 2100.0])
print("\n== FY2->M2 ==")
report(u2, M2, 270, 130, [4], [3.5], "单弹(t_l=4,τ=3.5)")
report(u2, M2, 270, 130, [4, 5, 6], [3.5, 3.5, 3.5], "三弹gap1同τ=3.5")
report(u2, M2, 270, 130, [4, 5], [3.5, 2.714], "两弹(τ=3.5/2.714) — MATLAB7.48s")
