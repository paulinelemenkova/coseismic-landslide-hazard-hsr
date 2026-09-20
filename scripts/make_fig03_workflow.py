from __future__ import annotations

import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
from matplotlib.font_manager import findfont, FontProperties

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})
_fp = FontProperties(family=plt.rcParams["font.sans-serif"])
assert "dejavu" not in findfont(_fp).lower(), "DejaVu fallback: install Nimbus Sans/Helvetica"

INK, MUTE = "#1f2328", "#3a4149"
GREEN  = dict(dark="#2E7D32", mid="#43A047", light="#EAF4EB")
BLUE   = dict(dark="#1565C0", mid="#1E88E5", light="#E8F1FB")
PURPLE = dict(dark="#6A1B9A", mid="#8E24AA", light="#F3EAF9")
ORANGE = dict(dark="#B45309", mid="#E08A1E", light="#FBEED9")
SLATE  = dict(dark="#37474F", mid="#546E7A", light="#ECEFF1")
TEAL   = dict(dark="#00695C", mid="#00897B", light="#E1F1EF")

XMIN, XMAX = 0.0, 136.0
LM, RM, GAP = 3.0, 3.0, 3.6
COLW = (XMAX - LM - RM - 3 * GAP) / 4.0
COLX = [LM + i * (COLW + GAP) for i in range(4)]
LH_BODY, LH_HEAD = 1.55, 1.8
PAD_T, PAD_B, GAP_HB, GAP_BLK, OUT_H = 0.9, 0.9, 0.55, 1.0, 6.0

def col_header_h(title: str) -> float:
    return max(6.2, 1.4 + (title.count("\n") + 1) * 2.0)

def block_h(lines: list[str]) -> float:
    return PAD_T + LH_HEAD + GAP_HB + len(lines) * LH_BODY + PAD_B

stages = [
 dict(c=GREEN, num="1", title="AI-BASED LANDSLIDE\nSUSCEPTIBILITY MODELLING", blocks=[
   ("INPUT \u2014 CONDITIONING FACTORS", ["Slope \u00b7 Lithology \u00b7 TWI", "Elevation \u00b7 Distance to fault"]),
   ("AI-ENHANCED MODELLING (AHP + ML)", ["AHP weighting \u2014 consistency check (CR = 0.072)", "AI assistance \u2014 feature importance", "& pattern learning"]),
   ("PLAUSIBILITY ASSESSMENT", ["Qualitative agreement with published", "inventories (G\u00f6r\u00fcm 2013; Duman 2005)"]),
 ], out="Landslide Susceptibility (LS)"),
 dict(c=BLUE, num="2", title="SEISMIC TRIGGERING\nASSESSMENT", blocks=[
   ("INPUT DATA", ["1999 \u0130zmit (Mw 7.6) scenario PGA field", "NAFZ fault geometry (Emre et al., 2013)"]),
   ("PROCESSING", ["PGA interpolation \u2192 distance-to-fault", "conditioning \u2192 Seismic Triggering (ST)", "CLH = LS \u00d7 ST \u2192 reclassification"]),
 ], out="Co-Seismic Landslide Hazard (CLH)"),
 dict(c=PURPLE, num="3", title="RAILWAY EXPOSURE\nANALYSIS", blocks=[
   ("INPUT DATA", ["Ankara\u2013\u0130stanbul HSR corridor (TCDD / OSM)", "Co-Seismic Landslide Hazard (CLH)"]),
   ("ANALYSIS", ["Buffer analysis (100 / 250 / 500 m) \u2192", "hazard overlay with corridor \u2192", "segment risk scoring & prioritization"]),
   ("OUTPUTS", ["Exposed sections \u00b7 high-risk segments", "\u00b7 risk ranking & statistics"]),
 ], out="Railway Exposure & Risk Map"),
 dict(c=ORANGE, num="4", title="CASCADING EARLY WARNING\nFRAMEWORK DESIGN", blocks=[
   ("1 \u00b7 Earthquake Early Warning (EEW)", ["P-wave detection & rapid magnitude"]),
   ("2 \u00b7 Ground-motion estimation", ["PGA / PGV nowcasting, study area"]),
   ("3 \u00b7 Dynamic landslide-hazard update", ["AI update of CLH from real-time", "PGA & fault proximity"]),
   ("4 \u00b7 Railway operational response", ["Speed reduction / stoppage / rerouting"]),
   ("5 \u00b7 Post-event inspection priority", ["High-risk segments inspected first"]),
 ], out="AI-Enhanced Cascading EW Framework"),
]
integration = ["Multi-source data integration", "AI-driven analytics & learning",
               "Real-time hazard visualization", "Automated alerts & notifications",
               "Historical database & feedback loop"]
outcomes = ["Improved safety of HSR infrastructure", "Risk-informed planning & maintenance",
            "Reduced response time via EEW", "Enhanced resilience in active-fault regions"]

def col_total(s: dict) -> float:
    hh = col_header_h(s["title"])
    bs = sum(block_h(l) for _, l in s["blocks"])
    return hh + bs + (len(s["blocks"]) + 1) * GAP_BLK + OUT_H + 1.0

STAGE_H = max(col_total(s) for s in stages)

Y_TITLE, Y_SUB1, Y_SUB2, Y_TOP = 98.0, 94.2, 91.6, 88.0
Y_BOT = Y_TOP - STAGE_H
MIDY = (Y_TOP + Y_BOT) / 2
BUS_Y = Y_BOT - 3.0
INT_TOP = Y_BOT - 5.5; INT_H = 8.6; INT_BOT = INT_TOP - INT_H
OUT_TOP = INT_BOT - 3.4; OUT_BH = 8.6; OUT_BOT = OUT_TOP - OUT_BH
CRED_Y = OUT_BOT - 3.0
YMIN, YMAX = CRED_Y - 1.8, 100.0
figw, figh = (XMAX - XMIN) * 0.1, (YMAX - YMIN) * 0.1

fig = plt.figure(figsize=(figw, figh))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(XMIN, XMAX); ax.set_ylim(YMIN, YMAX); ax.set_aspect("equal"); ax.axis("off")

def rbox(x0, y0, x1, y1, fc, ec, lw=1.0, rs=0.5, z=1):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
        boxstyle=f"round,pad=0,rounding_size={rs}", fc=fc, ec=ec, lw=lw, zorder=z, mutation_aspect=1.0))

def txt(x, y, s, size, color=INK, weight="normal", ha="left", va="center", style="normal", z=5):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va,
            style=style, zorder=z, linespacing=1.28)

txt((XMIN + XMAX) / 2, Y_TITLE, "METHODOLOGY WORKFLOW", 13.0, INK, "bold", ha="center")
txt((XMIN + XMAX) / 2, Y_SUB1,
    "Co-Seismic Landslide Hazard Assessment for High-Speed Railway Infrastructure in Northwestern T\u00fcrkiye:",
    9.2, MUTE, ha="center", style="italic")
txt((XMIN + XMAX) / 2, Y_SUB2,
    "A Knowledge-Driven AHP Model with an AI-Enhanced Cascading Early-Warning Framework",
    9.2, MUTE, ha="center", style="italic")

col_centers = []
for i, s in enumerate(stages):
    c = s["c"]; x0 = COLX[i]; x1 = x0 + COLW; inx0, inx1 = x0 + 1.1, x1 - 1.1
    col_centers.append((x0 + x1) / 2)
    rbox(x0, Y_BOT, x1, Y_TOP, c["light"], c["dark"], lw=1.3, rs=0.7, z=1)
    hh = col_header_h(s["title"]); htop, hbot = Y_TOP - 0.5, Y_TOP - hh + 0.4
    rbox(x0 + 0.6, hbot, x1 - 0.6, htop, c["dark"], c["dark"], lw=0, rs=0.5, z=2)
    cy = (htop + hbot) / 2
    ax.add_patch(Circle((x0 + 3.0, cy), 1.7, fc="white", ec=c["dark"], lw=1.0, zorder=3))
    txt(x0 + 3.0, cy, s["num"], 10.5, c["dark"], "bold", ha="center", z=4)
    txt(x0 + 5.6, cy, s["title"], 9.2, "white", "bold", z=4)

    out_inner_top = Y_BOT + 0.6 + OUT_H
    total_bh = sum(block_h(l) for _, l in s["blocks"]); n = len(s["blocks"])
    g = max((hbot - out_inner_top - total_bh) / (n + 1), 0.6)
    ycur = hbot - g
    for htitle, lines in s["blocks"]:
        bh = block_h(lines)
        rbox(inx0, ycur - bh, inx1, ycur, "white", c["mid"], lw=0.8, rs=0.4, z=3)
        ty = ycur - PAD_T
        txt(inx0 + 0.7, ty - LH_HEAD / 2, htitle, 8.4, c["dark"], "bold", z=5)
        by = ty - LH_HEAD - GAP_HB
        for ln in lines:
            txt(inx0 + 0.7, by - LH_BODY / 2, ln, 7.8, INK, z=5); by -= LH_BODY
        ycur -= bh + g
    rbox(x0 + 0.6, Y_BOT + 0.6, x1 - 0.6, Y_BOT + 0.6 + OUT_H, c["mid"], c["dark"], lw=0.8, rs=0.5, z=2)
    ofy = Y_BOT + 0.6 + OUT_H / 2
    txt((x0 + x1) / 2, ofy + 1.4, "OUTPUT", 7.6, "white", "bold", ha="center", z=4)
    txt((x0 + x1) / 2, ofy - 0.9, s["out"], 8.3, "white", "bold", ha="center", z=4)

for i in range(3):
    ax.add_patch(FancyArrowPatch((COLX[i] + COLW + 0.2, MIDY), (COLX[i + 1] - 0.2, MIDY),
        arrowstyle="-|>", mutation_scale=22, lw=3.0, color="#6b7280", zorder=2))

for xc in col_centers:
    ax.plot([xc, xc], [Y_BOT, BUS_Y], ls=(0, (4, 2.5)), lw=1.1, color="#9aa1a8", zorder=1)
ax.plot([col_centers[0], col_centers[-1]], [BUS_Y, BUS_Y], ls=(0, (4, 2.5)), lw=1.1, color="#9aa1a8", zorder=1)
xbus = (col_centers[0] + col_centers[-1]) / 2
ax.add_patch(FancyArrowPatch((xbus, BUS_Y), (xbus, INT_TOP + 0.1), arrowstyle="-|>",
    mutation_scale=18, lw=1.6, color="#6b7280", linestyle="dashed", zorder=1))

rbox(LM, INT_BOT, XMAX - RM, INT_TOP, SLATE["light"], SLATE["dark"], lw=1.3, rs=0.7, z=1)
txt((XMIN + XMAX) / 2, INT_TOP - 1.9,
    "INTEGRATION PLATFORM  \u00b7  AI-ENHANCED DECISION SUPPORT SYSTEM", 9.4, SLATE["dark"], "bold", ha="center")
row_y = INT_BOT + (INT_H - 3.0) / 2 + 0.2; seg = (XMAX - RM - LM) / len(integration)
for k, it in enumerate(integration):
    txt(LM + seg * (k + 0.5), row_y, it, 8.3, INK, ha="center")
    if k < len(integration) - 1:
        ax.plot([LM + seg * (k + 1)] * 2, [INT_BOT + 1.4, row_y + 2.0], lw=0.7, color=SLATE["mid"], zorder=2)
ax.add_patch(FancyArrowPatch(((XMIN + XMAX) / 2, INT_BOT - 0.1), ((XMIN + XMAX) / 2, OUT_TOP + 0.1),
    arrowstyle="-|>", mutation_scale=20, lw=2.4, color="#6b7280", zorder=2))

rbox(LM, OUT_BOT, XMAX - RM, OUT_TOP, TEAL["light"], TEAL["dark"], lw=1.3, rs=0.7, z=1)
txt((XMIN + XMAX) / 2, OUT_TOP - 1.9, "FINAL OUTCOMES", 9.4, TEAL["dark"], "bold", ha="center")
row_y2 = OUT_BOT + (OUT_BH - 3.0) / 2 + 0.2; seg2 = (XMAX - RM - LM) / len(outcomes)
for k, it in enumerate(outcomes):
    txt(LM + seg2 * (k + 0.5), row_y2, it, 8.3, INK, ha="center")
    if k < len(outcomes) - 1:
        ax.plot([LM + seg2 * (k + 1)] * 2, [OUT_BOT + 1.4, row_y2 + 2.0], lw=0.7, color=TEAL["mid"], zorder=2)

credit = ("Methodology workflow (knowledge-driven AHP; no landslide inventory or ROC/AUC step). "
          f"Software: Python {sys.version.split()[0]} \u00b7 Matplotlib {mpl.__version__}. "
          "Data: SRTM, GEM Global Active Faults, Natural Earth, TCDD/OSM. Source: authors.")
txt((XMIN + XMAX) / 2, CRED_Y, credit, 7.6, MUTE, ha="center", style="italic")

fig.savefig("fig03_workflow.pdf")
fig.savefig("fig03_workflow.png", dpi=600)
print("wrote fig03_workflow.pdf and fig03_workflow.png")
