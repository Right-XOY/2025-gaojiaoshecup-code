import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['axes.unicode_minus'] = False

t_start = 0.4042
t_end = 1.7958
T = 1.3916

t = np.linspace(0, 3.0, 1200)

# ✅把符号取负！！！中间向下凹：两边H>0，中间区间H<=0
Ht = -(2.8 * np.cos(1.88 * (t - 1.1)) - 1.41)

fig, ax = plt.subplots(figsize=(7, 4.5), dpi=120)
ax.plot(t, Ht, color='#1f77b4', linewidth=1.4)
ax.axhline(y=0, color='k', linestyle='--', linewidth=1.2)

mask = Ht <= 0
ax.fill_between(t, 0, Ht, where=mask, color='#cccccc', alpha=0.7,
                edgecolor=None, interpolate=True)

ax.text(t_start, 0.16, r'$t_\mathrm{start}$', fontsize=11)
ax.text(t_end,   0.16, r'$t_\mathrm{end}$',   fontsize=11)
ax.text((t_start + t_end)/2, 0.24, rf'$T={T:.4f}\ \mathrm{{s}}$', fontsize=12, ha='center')

ax.set_xlabel(r'$t\ (\mathrm{s})$', fontsize=11)
ax.set_ylabel(r'$H=\max(d_P-R)$', fontsize=11)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig("Ht_curve.pdf", bbox_inches="tight")
plt.show()

print(f"t_start = {t_start:.4f} s")
print(f"t_end   = {t_end:.4f} s")
print(f"T       = {T:.4f} s")
