import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
import sys
import os

out_prefix = sys.argv[1] if len(sys.argv) > 1 else "results/figures/shared_srna/srna_venn_pubstyle"
os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

# -----------------
# Counts
# -----------------
stationary_total = 519
exponential_total = 3752
shared = 24

stationary_only = stationary_total - shared   # 495
exponential_only = exponential_total - shared # 3728

# -----------------
# Style
# -----------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

fig = plt.figure(figsize=(8, 6))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_aspect("equal")
ax.axis("off")

# Canvas limits
ax.set_xlim(0.00, 1.28)
ax.set_ylim(0.00, 1.05)

# -----------------
# Colors
# -----------------
stationary_col   = (0.93, 0.76, 0.55, 0.85)  # peach
exponential_col  = (0.52, 0.65, 0.88, 0.85)  # blue
intersection_col = (0.40, 0.72, 0.40, 0.95)  # green

# -----------------
# Geometry
# -----------------
left_center  = (0.36, 0.62)
right_center = (0.83, 0.62)

left_r  = 0.22
right_r = 0.34

left_circle = Circle(
    left_center, left_r,
    facecolor=stationary_col,
    edgecolor="black",
    linewidth=1.2
)

right_circle = Circle(
    right_center, right_r,
    facecolor=exponential_col,
    edgecolor="black",
    linewidth=1.2
)

ax.add_patch(left_circle)
ax.add_patch(right_circle)

# True overlap by clipping left circle with right circle
intersection = Circle(
    left_center, left_r,
    facecolor=intersection_col,
    edgecolor="none"
)
intersection.set_clip_path(right_circle)
ax.add_patch(intersection)

# -----------------
# Numbers
# -----------------
# Left-only label
ax.text(
    left_center[0] - 0.11, left_center[1], f"{stationary_only}",
    ha="center", va="center",
    fontsize=26, fontweight="bold"
)

# Shared label: smaller and slightly left
shared_x = 0.54
ax.text(
    shared_x, left_center[1], f"{shared}",
    ha="center", va="center",
    fontsize=20, fontweight="bold"
)

# Right-only label
ax.text(
    right_center[0] + 0.15, right_center[1], f"{exponential_only}",
    ha="center", va="center",
    fontsize=26, fontweight="bold"
)

# -----------------
# Legend
# -----------------
handles = [
    Line2D([0], [0], marker='o', color='none',
           markerfacecolor=stationary_col, markeredgecolor='black',
           markeredgewidth=1.2, markersize=16,
           label="Stationary phase predicted sRNAs candidates"),
    Line2D([0], [0], marker='o', color='none',
           markerfacecolor=exponential_col, markeredgecolor='black',
           markeredgewidth=1.2, markersize=16,
           label="Exponential phase predicted sRNAs candidates"),
    Line2D([0], [0], marker='o', color='none',
           markerfacecolor=intersection_col, markeredgecolor='black',
           markeredgewidth=1.2, markersize=16,
           label="Shared loci"),
]

ax.legend(
    handles=handles,
    loc="lower left",
    bbox_to_anchor=(0.04, 0.03),
    frameon=False,
    fontsize=16,
    handletextpad=0.8
)

# -----------------
# Save
# -----------------
plt.savefig(out_prefix + ".png", dpi=600, bbox_inches="tight")
plt.savefig(out_prefix + ".pdf", bbox_inches="tight")
plt.close()
