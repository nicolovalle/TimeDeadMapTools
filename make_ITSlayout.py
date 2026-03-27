import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import argparse


def load_polygons(filename, onlyIB=False):
    polygons = []
    all_x = []
    all_y = []

    with open(filename, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()
            values = list(map(float, parts))

            bin_id = int(values[0])
            coords = values[1:]

            if onlyIB and bin_id >= 48:
                continue

            if len(coords) % 2 != 0:
                raise ValueError(
                    f"Line {line_num}: expected an even number of coordinates after the id, got {len(coords)}"
                )

            n_vertices = len(coords) // 2
            xs = coords[:n_vertices]
            ys = coords[n_vertices:]

            if onlyIB:
                xs = [i/3 for i in xs]
                ys = [i/3 for i in ys]

            pts = np.column_stack([xs, ys])

            polygons.append((bin_id, pts))
            all_x.extend(xs)
            all_y.extend(ys)

    return polygons, np.array(all_x), np.array(all_y)

def get_stave_id(binid):

    nstave = [0,12,28,48,72,102,144,192]
    for i in range(7):
        if nstave[i] <= binid < nstave[i+1]:
            return binid - nstave[i]

def plot_polygons(polygons, all_x, all_y, show_bin_labels=False, onlyIB=False):
    fig, ax = plt.subplots(figsize=(12, 12))

    scalelabel = 1.5 if onlyIB else 1

    for bin_id, pts in polygons:
        poly = Polygon(
            pts,
            closed=True,
            fill=False,
            edgecolor="black",
            linewidth=0.8
        )
        ax.add_patch(poly)

        if True: #show_bin_labels:
            center = pts.mean(axis=0)
            ax.text(
                center[0], center[1], f"{get_stave_id(bin_id)}\n({bin_id})",
                ha="center", va="center", fontsize=6*scalelabel
            )

    xmin, xmax = all_x.min(), all_x.max()
    ymin, ymax = all_y.min(), all_y.max()

    dx = xmax - xmin
    dy = ymax - ymin
    pad = 0.15 * max(dx, dy)

    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_aspect("equal")

    # Use the origin as the center of the polar reference
    cx, cy = 0.0, 0.0

    R = 1.08 * max(
        np.max(np.abs(all_x)),
        np.max(np.abs(all_y))
    )

    angles = [i*np.pi/16 for i in range(32)]
    """
        0,
        np.pi / 4,
        np.pi / 2,
        3 * np.pi / 4,
        np.pi,
        5 * np.pi / 4,
        3 * np.pi / 2,
        7 * np.pi / 4
    ]
"""
    labels = [f'{a:.2f}' for a in angles]
    """
        r"$0$",
        r"$\pi/4$",
        r"$\pi/2$",
        r"$3\pi/4$",
        r"$\pi$",
        r"$5\pi/4$",
        r"$3\pi/2$",
        r"$7\pi/4$"
    ]
"""

    theta = np.linspace(0, 2 * np.pi, 600)
    ax.plot(
        cx + R * np.cos(theta),
        cy + R * np.sin(theta),
        linestyle="--",
        linewidth=0.8,
        alpha=0.5
    )

    for ang, lab in zip(angles, labels):
        x_end = cx + R * np.cos(ang)
        y_end = cy + R * np.sin(ang)

        ax.plot(
            [cx, x_end],
            [cy, y_end],
            linestyle=":",
            linewidth=0.8,
            alpha=0.7
        )

        rt = 1.08 * R
        xt = cx + rt * np.cos(ang)
        yt = cy + rt * np.sin(ang)

        ax.text(
            xt, yt, lab,
            fontsize=11*scalelabel,
            ha="center", va="center"
        )

    ax.plot(cx, cy, "o", markersize=3)

    ax.set_title("Approximate position of ITS staves")
    slab = "" if onlyIB else "(IBx3)"
    ax.set_xlabel(f"x (mm) {slab}", fontsize = 12*scalelabel)
    ax.set_ylabel(f"y (mm) {slab}", fontsize = 12*scalelabel)
    ax.tick_params(axis='both', labelsize=12*scalelabel)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot THPoly-like polygons from txt file"
    )
    parser.add_argument("file", help="Input txt file")
    parser.add_argument("--labels", action="store_true", help="Show bin labels")
    parser.add_argument("--onlyib", action="store_true", help="Display only IB")
    args = parser.parse_args()

    polygons, all_x, all_y = load_polygons(args.file,onlyIB=args.onlyib)
    plot_polygons(polygons, all_x, all_y, show_bin_labels=args.labels, onlyIB=args.onlyib)
