import json
import matplotlib.pyplot as plt
from collections import defaultdict
import sys

if len(sys.argv[1:]) != 2:
    print(f'Usage: python3 {sys.argv[0]} <file1.json,file2.json> <spec>')
    sys.exit(1)

# Input files
infiles_str = str(sys.argv[1])
print(infiles_str)
if ',' in infiles_str:
    infiles = infiles_str.split(',')
else:
    infiles = [infiles_str]

spec = str(sys.argv[2])

clean_spec = ''.join(c for c in spec if c.isalpha() or c.isdigit())
print(f'Output will be saved in {clean_spec}.pdf')

filtered_data = []

# Read all files and tag each selected entry with its source file
for iff in infiles:
    with open(iff, "r") as f:
        print(f'Reading {iff}')
        data = json.load(f)

        filtered_data_iff = []
        for d in data:
            if d.get("bk_quality", "").lower() == "good" and not d.get("quality"):
                d_copy = d.copy()
                d_copy["_source_file"] = iff
                filtered_data_iff.append(d_copy)

        filtered_data.extend(filtered_data_iff)
        print(f'Selected runs so far: {len(filtered_data)}')

# Group by fill
runs_by_fill = defaultdict(list)
for d in filtered_data:
    fill = int(d["bk_fill"])
    runs_by_fill[fill].append(d)

# Sort fills
sorted_fills = sorted(runs_by_fill.keys())

# Determine the source file associated with each fill
# If a fill appears in more than one file, mark it as MIXED
fill_source = {}
for fill, runs in runs_by_fill.items():
    sources = {d["_source_file"] for d in runs}
    if len(sources) == 1:
        fill_source[fill] = next(iter(sources))
    else:
        fill_source[fill] = "MIXED"

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

# Plot run duration as filled area in background
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
ib_plot = ax1.plot(
    x_vals, ib_dtimes, 'o-',
    label="IB Dead Time",
    color='tab:blue',
    zorder=2
)
ob_plot = ax1.plot(
    x_vals, ob_dtimes, 's-',
    label="OB Dead Time",
    color='tab:orange',
    zorder=2
)

# Fill numbers on x-axis
xtick_positions = [x for x, _ in fill_labels]
xtick_labels = [label for _, label in fill_labels]
ax1.set_xticks(xtick_positions)
ax1.set_xticklabels(xtick_labels, rotation=90)

# Add vertical black lines at boundaries between different input files
for i in range(1, len(sorted_fills)):
    prev_fill = sorted_fills[i - 1]
    curr_fill = sorted_fills[i]

    if fill_source[prev_fill] != fill_source[curr_fill]:
        ax1.axvline(
            x=i - 0.5,
            color='red',
            linewidth=2,
            zorder=1
        )

# Labels and grid
ax1.set_xlabel("Fill Number")
ax1.set_ylabel("Dead Time [%]")
ax2.set_ylabel("Run Duration [min]")
ax1.grid(True)

# Combine legends from both axes
duration_proxy = ax2.fill_between([], [], color='gray', alpha=0.2, label='Run Duration')
lines = ib_plot + ob_plot + [duration_proxy]
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels, loc="upper left")

plt.title(f"IB/OB Dead Time per Fill - {spec}")
plt.tight_layout()
plt.savefig(f'{clean_spec}.pdf', dpi=300, bbox_inches='tight')
