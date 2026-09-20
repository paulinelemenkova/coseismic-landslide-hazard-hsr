import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Rectangle
import matplotlib.colors as mc
import matplotlib.font_manager as fm
import os

fp = "/usr/share/fonts/urw-base35/NimbusSans-Regular.otf"
fpb = "/usr/share/fonts/urw-base35/NimbusSans-Bold.otf"
if os.path.exists(fp):
    fm.fontManager.addfont(fp)
    if os.path.exists(fpb):
        fm.fontManager.addfont(fpb)
    plt.rcParams["font.family"] = fm.FontProperties(fname=fp).get_name()
plt.rcParams.update({"savefig.bbox": "tight"})

COL = ["#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7"]
INK, GREY = "#1a1a1a", "#555555"

def tint(hexc, a):
    r, g, b = mc.to_rgb(hexc)
    return (1 - a + a * r, 1 - a + a * g, 1 - a + a * b)

layers = [
    dict(t="Earthquake\nEarly Warning", d="Detect P-waves and\nissue an EEW alert.",
         i="\u2022 P-wave arrivals\n\u2022 Seismic network (KOERI)", p="EEW processing\n(real-time)",
         o="EEW alert:\nlocation, $M_w$,\nS-wave ETA, PGA range", lt="5\u201330 s", lt2="(epicentral distance)",
         k="Alert to railway\noperations control"),
    dict(t="Ground-Motion\nEstimation", d="Estimate PGA from EEW\nparameters and GMPEs.",
         i="\u2022 EEW magnitude & distance\n\u2022 Segment locations\n\u2022 Regional GMPEs", p="PGA estimation\n(rule-based)",
         o="Scenario / real-time\nPGA field (g)", lt="~2\u20135 s", lt2="after EEW",
         k="PGA values to\nhazard engine"),
    dict(t="Dynamic Landslide\nHazard Update", d="Compute CLH; flag\nthreshold exceedances.",
         i="\u2022 AHP susceptibility (LS)\n\u2022 Distance-to-fault\n\u2022 Real-time PGA (Layer 2)",
         p="CLH = LS \u00d7 ST\nST = min(PGA/PGA$_{ref}$, 1)", o="CLH per segment\n\u2192 5 hazard classes\n+ flagged segments",
         lt="~2\u20135 s", lt2="processing", k="List of flagged\nhigh-hazard segments"),
    dict(t="Railway Operational\nResponse", d="Automated actions on\nflagged segments.",
         i="\u2022 Flagged segments (Layer 3)\n\u2022 Train positions & speeds\n\u2022 Lead time; thresholds",
         p="Decision engine\n(operational logic)", o="Speed restriction /\nemergency braking /\nservice halt",
         lt="0\u201310 s", lt2="(within lead time)", k="Safety commands to\ntrains & wayside"),
    dict(t="Post-Event Inspection\nPrioritisation", d="Rank segments by CLH\nfor field inspection.",
         i="\u2022 Final CLH map\n\u2022 Observed ground motion\n\u2022 Infrastructure inventory",
         p="Prioritisation engine\n(ranking)", o="Inspection priority\n(Priority 1 \u2192 4)",
         lt="min\u2013hours", lt2="after event", k="Inspection plan to\nfield teams"),
]

fig = plt.figure(figsize=(13.2, 9.075))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 16); ax.set_ylim(0, 11); ax.set_aspect("equal"); ax.axis("off")

def rbox(x0, x1, y0, y1, fc, ec, lw=1.0, z=3, rs=0.06):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                 boxstyle=f"round,pad=0.01,rounding_size={rs}", fc=fc, ec=ec, lw=lw,
                 zorder=z, mutation_aspect=1.0))

def tx(x, y, s, fs, c=INK, w="normal", ha="center", va="center", z=6, ls=1.22, style="normal"):
    ax.text(x, y, s, fontsize=fs, color=c, fontweight=w, ha=ha, va=va, zorder=z,
            linespacing=ls, fontstyle=style)

def arrow(x0, y0, x1, y1, c, lw=1.8, z=5, ms=10, dash=None):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms,
                 lw=lw, color=c, zorder=z, shrinkA=0, shrinkB=0,
                 linestyle=dash if dash else "solid"))

tx(8.0, 10.62, "AI-Enhanced Cascading Early Warning Framework", 15.5, INK, "bold")
tx(8.0, 10.18, "for HSR Co-Seismic Landslide Risk Management (Istanbul\u2013Kocaeli corridor)", 9.5, GREY)
for xc, s in [(4.60, "INPUTS"), (7.20, "PROCESS"), (9.88, "OUTPUT"),
              (12.15, "LEAD TIME"), (14.45, "KEY OUTPUT / DECISION")]:
    tx(xc, 9.86, s, 8.2, GREY, "bold")

xL0, xL1 = 0.30, 3.22; xi0, xi1 = 3.42, 5.78; xp0, xp1 = 6.04, 8.36; xo0, xo1 = 8.62, 11.06
xt0, xt1 = 11.30, 13.00; xk0, xk1 = 13.18, 15.72; bh = 1.06
centers = [9.20, 7.80, 6.40, 5.00, 3.60]

for k, L in enumerate(layers):
    c = COL[k]; cy = centers[k]; y0 = cy - bh / 2; y1 = cy + bh / 2
    rbox(xL0, xk1, cy - 0.62, cy + 0.62, tint(c, 0.055), "none", z=1, rs=0.07)
    rbox(xL0, xL1, y0, y1, tint(c, 0.17), c, lw=1.1, z=3)
    ax.add_patch(Rectangle((xL0 + 0.02, y0 + 0.04), 0.09, bh - 0.08, fc=c, ec="none", zorder=4))
    ax.add_patch(Circle((0.72, cy), 0.235, fc=c, ec="white", lw=1.4, zorder=5))
    tx(0.72, cy, str(k + 1), 12, "white", "bold")
    tx(1.14, cy + 0.19, L["t"], 8.4, INK, "bold", ha="left", va="center", ls=1.12)
    tx(1.14, cy - 0.24, L["d"], 6.5, GREY, ha="left", va="center", ls=1.14)
    rbox(xi0, xi1, y0, y1, "white", c, lw=1.0, z=4); tx(xi0 + 0.12, cy, L["i"], 7.0, INK, ha="left", ls=1.28)
    rbox(xp0, xp1, y0, y1, tint(c, 0.10), c, lw=1.0, z=4); tx((xp0 + xp1) / 2, cy, L["p"], 7.4, INK, "bold", ls=1.3)
    rbox(xo0, xo1, y0, y1, "white", c, lw=1.0, z=4)
    tx((xo0 + xo1) / 2, cy + (0.12 if k == 2 else 0.0), L["o"], 7.0, INK, ls=1.28)
    arrow(xi1 + 0.02, cy, xp0 - 0.02, cy, c); arrow(xp1 + 0.02, cy, xo0 - 0.02, cy, c)
    rbox(xt0, xt1, y0, y1, tint(c, 0.10), c, lw=1.0, z=4)
    tx((xt0 + xt1) / 2, cy + 0.14, L["lt"], 9.6, c, "bold"); tx((xt0 + xt1) / 2, cy - 0.22, L["lt2"], 6.5, GREY)
    rbox(xk0, xk1, y0, y1, "white", c, lw=1.0, z=4); tx((xk0 + xk1) / 2, cy, L["k"], 7.2, INK, "bold", ls=1.28)

haz = ["#1a9850", "#91cf60", "#fee08b", "#fc8d59", "#d73027"]; cx = 8.78; cyy = centers[2] - 0.34
for j, hc in enumerate(haz):
    ax.add_patch(Rectangle((cx + j * 0.30, cyy), 0.26, 0.14, fc=hc, ec="0.35", lw=0.4, zorder=6))
tx(cx + len(haz) * 0.30 + 0.22, cyy + 0.07, "VL\u2192VH", 6.2, GREY, ha="left")

for k in range(4):
    arrow(0.72, centers[k] - 0.24, 0.72, centers[k + 1] + 0.24, GREY, lw=2.0, z=5, ms=11)

kb_y1, kb_y0 = 2.58, 1.00
rbox(xL0, xk1, kb_y0, kb_y1, tint("#7F7F7F", 0.10), "#5a5a5a", lw=1.1, z=3, rs=0.06)
tx((xL0 + xk1) / 2, kb_y1 - 0.19, "DATA & KNOWLEDGE BASE  (PRE-EVENT)", 9.4, INK, "bold")
tx((xL0 + xk1) / 2, kb_y1 - 0.40,
   "\u2191 feeds Layers 2\u20135      \u2193 refined by post-event outputs (weights, thresholds, SOPs)",
   6.8, GREY, style="italic")
kb = ["AHP susceptibility &\nconditioning factors", "Fault geometry\n(NAFZ)",
      "Historical earthquakes\n& ShakeMaps", "Published landslide\noccurrence (qualitative)",
      "HSR corridor &\ninfrastructure data", "Operational rules,\nthresholds & SOPs"]
n = len(kb); pad = 0.35; w = ((xk1 - xL0) - 2 * pad) / n; iy0, iy1 = kb_y0 + 0.14, kb_y1 - 0.58
for j, item in enumerate(kb):
    bx0 = xL0 + pad + j * w + 0.06; bx1 = xL0 + pad + (j + 1) * w - 0.06
    rbox(bx0, bx1, iy0, iy1, "white", "#8a8a8a", lw=0.8, z=4)
    ax.add_patch(Rectangle((bx0 + 0.10, (iy0 + iy1) / 2 - 0.20), 0.07, 0.40, fc=COL[j % 5], ec="none", zorder=5))
    tx((bx0 + bx1) / 2 + 0.10, (iy0 + iy1) / 2, item, 6.8, INK, ls=1.2)

for xf in (3.0, 6.5, 10.0):
    arrow(xf, kb_y1 + 0.02, xf, kb_y1 + 0.34, "#5a5a5a", lw=1.4, z=5, ms=9, dash=(0, (3, 2)))
arrow(15.30, centers[4] - 0.53, 15.30, kb_y1 + 0.02, COL[4], lw=1.5, z=5, ms=10, dash=(0, (4, 2)))

out = "/sandbox/output" if os.path.isdir("/sandbox/output") else "."
fig.savefig(f"{out}/fig07_ew_framework.pdf")
fig.savefig(f"{out}/fig07_ew_framework.png", dpi=600)
print("done")
