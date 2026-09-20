import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib.patches import Patch
import matplotlib.font_manager as fm
import os

fp = "/usr/share/fonts/urw-base35/NimbusSans-Regular.otf"
if os.path.exists(fp):
    fm.fontManager.addfont(fp)
    plt.rcParams["font.family"] = fm.FontProperties(fname=fp).get_name()
plt.rcParams.update({"font.size": 8.5, "axes.linewidth": 0.8,
                     "xtick.direction": "out", "ytick.direction": "out",
                     "savefig.bbox": "tight"})
OI = {"orange": "#E69F00", "green": "#009E73", "blue": "#0072B2",
      "verm": "#D55E00", "grey": "#7F7F7F"}

factors = ["Slope angle", "Distance to fault", "Lithology", "TWI", "Elevation", "Land cover"]
w = np.array([0.312, 0.241, 0.178, 0.112, 0.094, 0.063])
assert abs(w.sum() - 1.0) < 1e-9

rng = np.random.default_rng(20250919)
N = 20000
wp = w * (1.0 + rng.uniform(-0.20, 0.20, size=(N, w.size)))
wp /= wp.sum(axis=1, keepdims=True)
p5, p95 = np.percentile(wp, [5, 95], axis=0)
base_order = np.argsort(-w)
top2 = set(base_order[:2].tolist())
p_top2 = np.mean([set(np.argsort(-r)[:2].tolist()) == top2 for r in wp])
p_full = np.mean([np.array_equal(np.argsort(-r), base_order) for r in wp])

thr = np.full(w.size, np.nan)
for pos, idx in enumerate(base_order):
    c = []
    if pos > 0:
        c.append(abs(w[base_order[pos - 1]] / w[idx] - 1))
    if pos < w.size - 1:
        c.append(abs(w[base_order[pos + 1]] / w[idx] - 1))
    thr[idx] = min(c) * 100

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.4, 3.7))
y = np.arange(w.size)[::-1]
tier_col = [OI["verm"], OI["verm"], OI["blue"], OI["blue"], OI["grey"], OI["grey"]]

axA.barh(y, w, height=0.62, color=tier_col, edgecolor="black", linewidth=0.5, zorder=3)
axA.errorbar(w, y, xerr=[w - p5, p95 - w], fmt="none", ecolor="black",
             elinewidth=1.2, capsize=3, capthick=1.2, zorder=4)
for yi, wi, hi in zip(y, w, p95):
    axA.text(hi + 0.006, yi, f"{wi:.3f}", va="center", ha="left", fontsize=8, zorder=5)
axA.set_yticks(y); axA.set_yticklabels(factors)
axA.set_xlim(0, 0.40); axA.set_xlabel("Normalised AHP weight")
axA.xaxis.set_major_locator(MultipleLocator(0.1)); axA.xaxis.set_minor_locator(AutoMinorLocator(2))
axA.yaxis.set_minor_locator(AutoMinorLocator(1))
axA.grid(axis="x", which="major", color="0.8", linewidth=0.6, zorder=0); axA.set_axisbelow(True)
axA.set_title("Factor weights with \u00b120% robustness envelope", fontsize=8.5, pad=4)
leg = axA.legend(handles=[Patch(fc=OI["verm"], ec="black", lw=.5, label="Dominant (rank 1\u20132)"),
                          Patch(fc=OI["blue"], ec="black", lw=.5, label="Moderate (rank 3\u20134)"),
                          Patch(fc=OI["grey"], ec="black", lw=.5, label="Low (rank 5\u20136)")],
                 loc="lower right", fontsize=6.6, frameon=True, framealpha=0.95,
                 handlelength=1.1, borderpad=0.5)
leg.get_frame().set_edgecolor("0.7"); leg.get_frame().set_linewidth(0.5)
axA.text(0.395, 2.15, f"Monte-Carlo (N={N:,}, \u00b120%):\nP(top-2 order kept) = {p_top2:.2f}\n"
                      f"P(full order kept) = {p_full:.2f}", ha="right", va="center", fontsize=6.8,
         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", lw=0.5), zorder=6)

thr_plot = thr[base_order]
fac_plot = [factors[i] for i in base_order]
sens = thr_plot < 20.0
col_plot = [OI["orange"] if s else OI["green"] for s in sens]
axB.barh(y, thr_plot, height=0.62, color=col_plot, edgecolor="black", linewidth=0.5, zorder=3)
for yi, ti, s in zip(y, thr_plot, sens):
    if s:
        axB.text(ti - 0.9, yi, f"{ti:.0f}%", va="center", ha="right", fontsize=8,
                 color="white", fontweight="bold", zorder=5)
    else:
        axB.text(ti + 0.8, yi, f"{ti:.0f}%", va="center", ha="left", fontsize=8, zorder=5)
axB.axvline(20, color="black", lw=1.2, ls=(0, (4, 2)), zorder=2)
axB.set_yticks(y); axB.set_yticklabels(fac_plot)
axB.set_xlim(0, 58); axB.set_xlabel("Min. weight change to alter rank (%)")
axB.xaxis.set_major_locator(MultipleLocator(20)); axB.xaxis.set_minor_locator(AutoMinorLocator(2))
axB.yaxis.set_minor_locator(AutoMinorLocator(1))
axB.grid(axis="x", which="major", color="0.8", linewidth=0.6, zorder=0); axB.set_axisbelow(True)
axB.set_title("Rank-reversal threshold per factor", fontsize=8.5, pad=4)
legB = axB.legend(handles=[Patch(fc=OI["green"], ec="black", lw=.5, label="Robust ($\\geq$20%)"),
                           Patch(fc=OI["orange"], ec="black", lw=.5, label="Sensitive ($<$20%)")],
                  loc="center right", bbox_to_anchor=(0.995, 0.40), fontsize=6.6, frameon=True,
                  framealpha=0.97, handlelength=1.1, borderpad=0.5)
legB.get_frame().set_edgecolor("0.7"); legB.get_frame().set_linewidth(0.5)

for ax, tag in ((axA, "(a)"), (axB, "(b)")):
    ax.text(0.985, 0.985, tag, transform=ax.transAxes, ha="right", va="top",
            fontsize=9, fontweight="bold")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(which="both", length=3.5, width=0.8)
    ax.tick_params(which="minor", length=2)

fig.tight_layout(w_pad=1.6)
out = "/sandbox/output" if os.path.isdir("/sandbox/output") else "."
fig.savefig(f"{out}/fig04_ahp_sensitivity.pdf")
fig.savefig(f"{out}/fig04_ahp_sensitivity.png", dpi=600)
print("thresholds %:", dict(zip(factors, np.round(thr, 1))))
print(f"P(top-2 kept)={p_top2:.3f}  P(full order kept)={p_full:.3f}")
