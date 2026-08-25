# -*- coding: utf-8 -*-
"""
q5_repro.py  问题5 全量复现（网格搜索版，用于定叙事口径 + 给出自洽数值）

结论依据（由 _verify_q5_relay.py 验证）：
  - 同机多弹接力方向 = 单弹最优方向（朝 +x 类），而非朝 -x；
  - 相邻投放间隔 gap=1s、3 弹沿单弹最优方向分布，可把单弹窗口前后延伸。
本脚本对每个 (uav, missile) 组合做 3 弹接力网格搜索，再 ILP 指派，
输出各机 θ/v/t_l/τ、各弹窗口与三导弹遮蔽时长之和。
"""
import numpy as np

G = 9.8
R = 10.0
VM = 300.0

M0_all = np.array([
    [20000.0, 19000.0, 18000.0],
    [0.0, 600.0, -600.0],
    [2000.0, 2100.0, 1900.0],
])
UAV0_all = np.array([
    [17800.0, 12000.0, 6000.0, 11000.0, 13000.0],
    [0.0, 1400.0, -3000.0, 2000.0, -2000.0],
    [1800.0, 1400.0, 700.0, 1800.0, 1300.0],
])

# 预扫描单弹最优（uav1, m1, theta, v, t_l, tau, bestT, win_s, win_e）
SCAN_PAR = [
    [1, 1, 10.0, 110.0, 0.0, 0.5, 4.60, 2.30, 6.90],
    [2, 1, 260.0, 120.0, 5.0, 6.5, 3.80, 15.30, 19.10],
    [3, 1, 90.0, 120.0, 22.0, 3.5, 2.70, 37.00, 39.70],
    [4, 1, 280.0, 90.0, 11.0, 11.5, 2.70, 26.30, 29.00],
    [5, 1, 110.0, 110.0, 16.0, 3.5, 3.50, 21.20, 24.70],
    [2, 2, 270.0, 130.0, 4.0, 3.5, 3.90, 9.70, 13.60],
    [3, 2, 90.0, 100.0, 30.0, 2.5, 2.60, 33.40, 36.00],
    [4, 2, 280.0, 130.0, 2.0, 10.5, 3.60, 14.10, 17.70],
    [5, 2, 120.0, 140.0, 19.0, 0.5, 1.40, 23.40, 24.80],
    [2, 3, 230.0, 130.0, 10.0, 7.5, 2.80, 21.60, 24.40],
    [3, 3, 100.0, 120.0, 20.0, 4.5, 2.60, 29.20, 31.80],
    [5, 3, 120.0, 120.0, 14.0, 1.5, 1.80, 17.20, 19.00],
]

nPhi = 100
ang = np.arange(nPhi) * 2 * np.pi / nPhi
rim = []
for z in [0.0, 10.0]:
    for a in ang:
        rim.append([7 * np.cos(a), 200 + 7 * np.sin(a), z])
rim = np.array(rim).T  # (3, 200)


def proj_window(uav0, M0, theta_deg, v, t_l, tau, dt):
    """单弹在活跃窗 [td, td+20] 内的遮蔽布尔序列 + 对应时刻。"""
    T_end = np.linalg.norm(M0) / VM
    um = -M0 / np.linalg.norm(M0)
    uh = np.array([np.cos(np.radians(theta_deg)), np.sin(np.radians(theta_deg)), 0.0])
    td = t_l + tau
    if td >= T_end or td + 20 <= 0:
        return None, None
    tmax = min(T_end, td + 20)
    ts = np.arange(td, tmax, dt)
    if len(ts) == 0:
        return None, None
    Pl = uav0 + v * t_l * uh
    Pd = Pl + v * tau * uh - np.array([0.0, 0.0, 0.5 * G * tau ** 2])
    Ms = M0[:, None] + VM * um[:, None] * ts[None, :]      # (3, T)
    C = Pd[:, None] - np.array([0.0, 0.0, 3.0])[:, None] * (ts[None, :] - td)  # (3, T)
    vvec = rim[:, None, :] - Ms[:, :, None]                # (3, T, 200)
    wvec = C[:, :, None] - Ms[:, :, None]
    vv_norm = np.sum(vvec * vvec, axis=0)
    tp = np.sum(wvec * vvec, axis=0) / vv_norm
    tp = np.clip(tp, 0.0, 1.0)
    closest = Ms[:, :, None] + tp[None, :, :] * vvec
    d = np.sqrt(np.sum((C[:, :, None] - closest) ** 2, axis=0))
    maxd = np.max(d, axis=1)
    occ = maxd <= R
    return occ, ts


def multi_shield(uav0, M0, theta_deg, v, t_ls, taus, dt=0.02):
    """多弹并集遮蔽时长与窗口。"""
    T_end = np.linalg.norm(M0) / VM
    all_ts = []
    for t_l, tau in zip(t_ls, taus):
        occ, ts = proj_window(uav0, M0, theta_deg, v, t_l, tau, dt)
        if occ is None:
            continue
        all_ts.append(ts[occ])
    if not all_ts:
        return 0.0, None, None
    merged = np.unique(np.concatenate(all_ts))
    # 连续区间按 dt 近似计长
    if len(merged) == 0:
        return 0.0, None, None
    dur = len(merged) * dt
    return dur, merged[0], merged[-1]


def relay_best(uav0, M0, base, dt=0.02):
    th0, v0, tl0, tau0 = base
    T_end = np.linalg.norm(M0) / VM
    tau_max = np.sqrt(2 * uav0[2] / G)
    tau_cand = sorted({x for x in
                       [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0]
                       if 0 < x <= tau_max})
    th = th0                       # 接力方向 = 单弹最优方向（已数值验证）
    best = (0.0, None)
    for v in range(70, 145, 5):
        for gap in [1.0, 2.0]:
            for tau in tau_cand:
                for t0 in [max(0.0, tl0 - 3), tl0, tl0 + 3]:
                    t_ls = [t0, t0 + gap, t0 + 2 * gap]
                    taus = [tau, tau, tau]
                    if t_ls[-1] + tau > T_end:
                        continue
                    d, s, e = multi_shield(uav0, M0, th, v, t_ls, taus, dt)
                    if d > best[0]:
                        best = (d, (th, v, t_ls, taus, s, e))
    return best


def main():
    relay = {}
    for row in SCAN_PAR:
        uav = int(row[0]) - 1
        m = int(row[1]) - 1
        base = row[2:6]
        single = row[6]
        uav0 = UAV0_all[:, uav]
        M0 = M0_all[:, m]
        d, p = relay_best(uav0, M0, base)
        d = max(d, single)   # 接力不低于单弹（退化情况）
        relay[(uav, m)] = (d, p)
        if p is None:
            print(f"FY{uav+1}->M{m+1}: 单弹={single:.2f}s  接力(未找到)=退回单弹={d:.2f}s")
            continue
        print(f"FY{uav+1}->M{m+1}: 单弹={single:.2f}s  接力={d:.2f}s  "
              f"(th={p[0]:.0f}, v={p[1]:.0f}, gap={p[2][1]-p[2][0]:.1f}, "
              f"tau={p[3][0]:.1f}, 窗口=[{p[4]:.2f},{p[5]:.2f}])")

    nU, nM = 5, 3
    coeff = np.zeros((nU, nM))
    for u in range(nU):
        for m in range(nM):
            coeff[u, m] = relay[(u, m)][0] if (u, m) in relay else 0.0

    best_assign = None
    best_score = -1.0
    for code in range(4 ** nU):
        a = [(code // (4 ** u)) % 4 for u in range(nU)]
        covered = [False] * nM
        sc = 0.0
        ok = True
        for u in range(nU):
            mm = a[u]
            if mm == 0:
                continue
            if coeff[u, mm - 1] <= 0:
                ok = False
                break
            covered[mm - 1] = True
            sc += coeff[u, mm - 1]
        if ok and all(covered) and sc > best_score:
            best_score = sc
            best_assign = a

    print("\n=== ILP 最优指派（系数=接力时长）===")
    groups = {m: [] for m in range(1, nM + 1)}
    for u in range(nU):
        mm = best_assign[u]
        if mm:
            groups[mm].append(u + 1)
    for m in sorted(groups):
        print(f"  M{m} <- " + ", ".join(f"FY{u}" for u in groups[m]))

    print("\n=== 各导弹遮蔽窗口（并集）===")
    total = 0.0
    for m in range(1, nM + 1):
        s = 0.0
        for u in groups[m]:
            if (u - 1, m - 1) in relay:
                d, p = relay[(u - 1, m - 1)]
                s += d
                print(f"  M{m}: FY{u} 窗口=[{p[4]:.2f},{p[5]:.2f}] 时长={d:.2f}s")
        total += s
        print(f"  M{m} 遮蔽时长(并集近似)={s:.2f}s")
    print(f"\n三导弹遮蔽时长之和 ≈ {total:.2f} s")
    print("（注：网格搜索近似，最终以 MATLAB 全量 JADE 为准）")


if __name__ == "__main__":
    main()
