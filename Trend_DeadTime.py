import json
import matplotlib.pyplot as plt
from collections import defaultdict
import sys

if len(sys.argv[1:]) != 2:
    print(f'Usage: python3 {sys.argv[0]} <file.json> <spec>')
    exit()
# Load JSON

infile = str(sys.argv[1])
spec = str(sys.argv[2])

with open(infile, "r") as f:
    data = json.load(f)

# Filter for good quality runs
filtered_data = [
    d for d in data
    if d.get("bk_quality", "").lower() == "good" and not d.get("quality")
]

# Group by fill
runs_by_fill = defaultdict(list)
for d in filtered_data:
    fill = int(d["bk_fill"])
    runs_by_fill[fill].append(d)

# Sort fills
sorted_fills = sorted(runs_by_fill.keys())

# Prepare plotting data
x_vals = []
ib_dtimes = []
ob_dtimes = []
rct_durations = []
fill_labels = []

for fill_idx, fill in enumerate(sorted_fills):
    runs = runs_by_fill[fill]
    for offset, d in enumerate(runs):
        x = fill_idx + offset * 0.15
        x_vals.append(x)
        ib_dtimes.append(float(d["ib_dtime"]))
        ob_dtimes.append(float(d["ob_dtime"]))
        rct_durations.append(float(d["rct_duration"]))
    fill_labels.append((fill_idx + (len(runs) - 1) * 0.075, str(fill)))

# Start plotting
fig, ax1 = plt.subplots(figsize=(24, 6))

# Secondary axis for duration
ax2 = ax1.twinx()

# Plot run duration as filled area (shadow), in background
ax2.fill_between(
    x_vals,
    rct_durations,
    color='gray',
    alpha=0.2,
    step='mid',
    label='Run Duration',
    zorder=0,
)

# Plot IB and OB
ib_plot = ax1.plot(x_vals, ib_dtimes, 'o-', label="IB Dead Time", color='tab:blue', zorder=2)
ob_plot = ax1.plot(x_vals, ob_dtimes, 's-', label="OB Dead Time", color='tab:orange', zorder=2)

# Fill numbers on x-axis
xtick_positions = [x for x, _ in fill_labels]
xtick_labels = [label for _, label in fill_labels]
ax1.set_xticks(xtick_positions)
ax1.set_xticklabels(xtick_labels, rotation=90)

# Labels and grid
ax1.set_xlabel("Fill Number")
ax1.set_ylabel("Dead Time [%]")
ax2.set_ylabel("Run Duration [min]")
ax1.grid(True)

# Combine legends from both axes
lines = ib_plot + ob_plot + [ax2.fill_between([], [], color='gray', alpha=0.2, label='Run Duration')]
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels, loc="upper left")

plt.title(f"IB/OB Dead Time per Fill - {spec}")
plt.tight_layout()
#plt.show()
plt.savefig(infile.replace('.json','.pdf'), dpi=300, bbox_inches='tight')
