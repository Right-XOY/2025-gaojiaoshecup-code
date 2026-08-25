import numpy as np

G = 9.8
R = 10.0
VM = 300.0
M0 = np.array([20000.0, 0.0, 2000.0])
um = -M0 / np.linalg.norm(M0)
T_end = np.linalg.norm(M0) / VM

nPhi = 100
ang = np.arange(nPhi) * 2 * np.pi / nPhi
rim = []
for z in [0.0, 10.0]:
    for a in ang:
        rim.append([7 * np.cos(a), 200 + 7 * np.sin(a), z])
rim = np.array(rim).T  # (3, 200)

def shield(uav0, theta_deg, v, t_ls, taus, dt=0.005):
    uh = np.array([np.cos(np.radians(theta_deg)), np.sin(np.radians(theta_deg)), 0.0])
    ts = np.arange(0.0, T_end, dt)
    N = len(ts)
    Ms = M0[:, None] + VM * um[:, None] * ts[None, :]  # (3, N)
    shielded = np.zeros(N, dtype=bool)
    for t_l, tau in zip(t_ls, taus):
        Pl = uav0 + v * t_l * uh
        Pd = Pl + v * tau * uh - np.array([0.0, 0.0, 0.5 * G * tau**2])
        td = t_l + tau
        # 云团中心 (3, N)
        sink = np.array([0.0, 0.0, 3.0])[:, None] * (ts[None, :] - td)
        C = Pd[:, None] - sink
        # 距离：C(3,N) 到 线段 M(3,N)-P(3,200)，统一用 (3,N,200) 布局
        vvec = rim[:, None, :] - Ms[:, :, None]            # (3,N,200)
        wvec = C[:, :, None] - Ms[:, :, None]              # (3,N,200)
        vv_norm = np.sum(vvec * vvec, axis=0)              # (N,200)
        tp = np.sum(wvec * vvec, axis=0) / vv_norm         # (N,200)
        tp = np.clip(tp, 0.0, 1.0)
        closest = Ms[:, :, None] + tp[None, :, :] * vvec   # (3,N,200)
        d = np.sqrt(np.sum((C[:, :, None] - closest) ** 2, axis=0))  # (N,200)
        maxd = np.max(d, axis=1)                            # (N,)
        active = (ts >= td) & (ts <= td + 20)
        shielded |= active & (maxd <= R)
    dur = np.sum(shielded) * dt
    if shielded.any():
        idx = np.where(shielded)[0]
        return dur, ts[idx[0]], ts[idx[-1]]
    return dur, None, None

uav0 = np.array([17800.0, 0.0, 1800.0])

print("=== 验证判据：单弹最优 FY1->M1 theta=10 v=110 t_l=0 tau=0.5 ===")
d, s, e = shield(uav0, 10.0, 110.0, [0.0], [0.5])
print(f"  时长={d:.4f}s  窗口=[{s:.2f},{e:.2f}]  (预扫描应为 4.60, [2.30,6.90])")

print("=== 朝-x 单弹 theta=180 v=110 t_l=0 tau=0.5 ===")
d, s, e = shield(uav0, 180.0, 110.0, [0.0], [0.5])
print(f"  时长={d:.4f}s  窗口=[{s:.2f},{e:.2f}]")

print("=== 朝-x 三弹接力 theta=180 v=110 t_l=[0,1,2] tau=[0.5]*3 ===")
d, s, e = shield(uav0, 180.0, 110.0, [0.0, 1.0, 2.0], [0.5, 0.5, 0.5])
print(f"  时长={d:.4f}s  窗口=[{s:.2f},{e:.2f}]")

print("=== 扫描朝-x 接力：v 与 t_l 间隔 ===")
best = (0, None)
for v in [100.0, 110.0, 120.0, 130.0, 140.0]:
    for gap in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]:
        d, s, e = shield(uav0, 180.0, v, [0.0, gap, 2 * gap], [0.5, 0.5, 0.5])
        if d > best[0]:
            best = (d, (v, gap, s, e))
print(f"  朝-x 三弹最佳：时长={best[0]:.4f}s  (v={best[1][0]}, gap={best[1][1]}, 窗口=[{best[1][2]:.2f},{best[1][3]:.2f}])")

print("=== 扫描朝+x 接力：v 与 t_l 间隔（弹沿+x，遮蔽前移）===")
best2 = (0, None)
for v in [90.0, 100.0, 110.0, 120.0]:
    for gap in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]:
        d, s, e = shield(uav0, 10.0, v, [0.0, gap, 2 * gap], [0.5, 0.5, 0.5])
        if d > best2[0]:
            best2 = (d, (v, gap, s, e))
print(f"  朝+x 三弹最佳：时长={best2[0]:.4f}s  (v={best2[1][0]}, gap={best2[1][1]}, 窗口=[{best2[1][2]:.2f},{best2[1][3]:.2f}])")

print("=== 朝+x 三弹接力（合法 gap>=1）：v × gap × τ 联合扫描 ===")
best5 = (0, None)
for v in [85.0, 90.0, 95.0, 100.0, 105.0, 110.0]:
    for gap in [1.0, 1.5, 2.0, 2.5, 3.0]:
        for tau in [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
            d, s, e = shield(uav0, 10.0, v, [0.0, gap, 2 * gap], [tau, tau, tau])
            if d > best5[0]:
                best5 = (d, (v, gap, tau, s, e))
print(f"  合法朝+x 三弹最佳：时长={best5[0]:.4f}s  (v={best5[1][0]}, gap={best5[1][1]}, tau={best5[1][2]}, 窗口=[{best5[1][3]:.2f},{best5[1][4]:.2f}])")
print(f"  => 相比单弹 4.505s 提升 {best5[0]-4.505:.4f}s ({(best5[0]-4.505)/4.505*100:.1f}%)")



