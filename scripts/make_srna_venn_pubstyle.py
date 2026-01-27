import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

# -----------------
# Counts
# -----------------
stationary_total = 964
exponential_total = 8742
shared = 104

stationary_only = stationary_total - shared
exponential_only = exponential_total - shared

# -----------------
# Style
# -----------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

fig = plt.figure(figsize=(6, 4.8))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_aspect("equal")
ax.axis("off")

# -----------------
# Colors (match your current figure)
# -----------------
stationary_col   = (0.98, 0.78, 0.55, 0.85)  # orange
exponential_col  = (0.52, 0.65, 0.88, 0.85)  # blue
intersection_col = (0.40, 0.72, 0.40, 0.95)  # green

# -----------------
# Geometry
# -----------------
left_center  = (0.40, 0.55)
right_center = (0.67, 0.55)
left_r  = 0.24
right_r = 0.32

left_circle = Circle(left_center, left_r,
                     facecolor=stationary_col,
                     edgecolor="black", linewidth=1.0)
right_circle = Circle(right_center, right_r,
                      facecolor=exponential_col,
                      edgecolor="black", linewidth=1.0)

ax.add_patch(left_circle)
ax.add_patch(right_circle)

# TRUE intersection (clip the left circle by the right circle)
intersection = Circle(left_center, left_r,
                      facecolor=intersection_col,
                      edgecolor="none")
intersection.set_clip_path(right_circle)
ax.add_patch(intersection)

# -----------------
# Numbers
# -----------------
ax.text(left_center[0] - 0.12, left_center[1], f"{stationary_only}",
        ha="center", va="center", fontsize=15, fontweight="bold")
ax.text((left_center[0] + right_center[0]) / 2, left_center[1], f"{shared}",
        ha="center", va="center", fontsize=15, fontweight="bold")
ax.text(right_center[0] + 0.17, right_center[1], f"{exponential_only}",
        ha="center", va="center", fontsize=15, fontweight="bold")

# -----------------
# Legend (no title)
# -----------------
handles = [
    Line2D([0], [0], marker='o', color='none',
           markerfacecolor=stationary_col, markeredgecolor='black',
           markersize=14, label="Stationary phase sRNAs"),
    Line2D([0], [0], marker='o', color='none',
           markerfacecolor=exponential_col, markeredgecolor='black',
           markersize=14, label="Exponential phase sRNAs"),
    Line2D([0], [0], marker='o', color='none',
           markerfacecolor=intersection_col, markeredgecolor='black',
           markersize=14, label="Shared sRNAs"),
]
ax.legend(handles=handles, loc="lower left",
          bbox_to_anchor=(0.03, 0.03),
          frameon=False)

# -----------------
# Save
# -----------------
plt.savefig("srna_venn_no_title.png", dpi=600, bbox_inches="tight")
plt.savefig("srna_venn_no_title.pdf", bbox_inches="tight")
plt.show()
