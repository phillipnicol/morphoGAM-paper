#!/usr/bin/env python

"""Draft Spateo digitization workflow for the granule simulation geometry.

Run from analysis/granule_t_sim/code, for example:

    conda run -n spateo-release-morphogam \
      python spateo_granule_digitize.py

The script mirrors analysis/ca3_newsim/code/spateo_ca3_digitize*.ipynb:
it loads granule_xy.csv, creates PCA-aligned image coordinates, builds a
contour, picks four Spateo corner points, runs st.dd.digitize(), and writes
granule_xy_digitized.csv for run_spateo_models.py.
"""

import argparse
import os
import warnings
from pathlib import Path

# Set writable cache locations before importing Spateo/Numba-backed modules.
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/cache-spateo")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mpl-spateo")
os.environ.setdefault("NUMBA_CACHE_DIR", "/private/tmp/numba-spateo")

import anndata as ad
import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spateo as st

# Spateo's digitize() currently emits many pandas chained-assignment
# FutureWarnings while filling adata.obs. They are noisy but not fatal here.
warnings.filterwarnings("ignore", category=FutureWarning, module="spateo.digitization.grid")


CORNER_NAMES = ("pnt_xy", "pnt_xY", "pnt_Xy", "pnt_XY")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Digitize granule coordinates with Spateo's layer/column workflow."
    )
    parser.add_argument("--data-dir", type=Path, default=Path("../data"))
    parser.add_argument("--plot-dir", type=Path, default=Path("../plots"))
    parser.add_argument("--xy-file", type=Path, default=None)
    parser.add_argument("--output-file", type=Path, default=None)
    parser.add_argument("--plot-prefix", default="granule_spateo_digitize")
    parser.add_argument("--digitize-pad", type=int, default=50)
    parser.add_argument("--digitize-scale", type=int, default=120)
    parser.add_argument("--max-itr", type=int, default=5000)
    parser.add_argument("--lh", type=float, default=1)
    parser.add_argument("--hh", type=float, default=100)
    parser.add_argument(
        "--flip-longitudinal",
        action="store_true",
        help="Reverse the PCA longitudinal/image-column direction.",
    )
    parser.add_argument(
        "--flip-radial",
        action="store_true",
        help="Reverse the PCA radial/image-layer direction.",
    )
    parser.add_argument(
        "--corner-rows",
        nargs="+",
        default=None,
        metavar="NAME=ROW",
        help=(
            "Optional zero-based row overrides for corner spots, e.g. "
            "pnt_xy=10 pnt_xY=42 pnt_Xy=77 pnt_XY=123. "
            "Rows are snapped to the nearest ROI contour point."
        ),
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Do not write corner/digitized coordinate QC plots.",
    )
    return parser.parse_args()


def normalize01(x):
    x = np.asarray(x, dtype=float)
    rng = x.max() - x.min()
    if rng == 0:
        return np.zeros_like(x, dtype=float)
    return (x - x.min()) / rng


def parse_corner_rows(values):
    if values is None:
        return None

    parsed = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"Corner row must have NAME=ROW form, got {value!r}")
        name, row = value.split("=", 1)
        if name not in CORNER_NAMES:
            raise ValueError(f"Unknown corner name {name!r}; expected one of {CORNER_NAMES}")
        parsed[name] = int(row)

    missing = set(CORNER_NAMES) - set(parsed)
    if missing:
        raise ValueError(f"Missing corner row override(s): {sorted(missing)}")
    return parsed


def load_coordinates(xy_file):
    xy = pd.read_csv(xy_file)
    required_cols = ["xcoord", "ycoord"]
    missing_cols = [col for col in required_cols if col not in xy.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in {xy_file}: {missing_cols}")

    coords = xy[required_cols].to_numpy(dtype=float)
    adata = ad.AnnData(X=np.ones((coords.shape[0], 1)))
    adata.obs_names = [f"spot_{i + 1}" for i in range(coords.shape[0])]
    adata.var_names = ["dummy"]
    adata.obsm["spatial"] = coords
    adata.obs["xcoord"] = coords[:, 0]
    adata.obs["ycoord"] = coords[:, 1]
    adata.uns["__type"] = "UMI"
    return xy, coords, adata


def pca_image_coordinates(coords, args):
    center = coords.mean(axis=0)
    coords_centered = coords - center
    _, _, vt = np.linalg.svd(coords_centered, full_matrices=False)

    longitudinal_axis = vt[0]
    radial_axis = vt[1]
    longitudinal_score = coords_centered @ longitudinal_axis
    radial_score = coords_centered @ radial_axis

    if args.flip_longitudinal:
        longitudinal_score = -longitudinal_score
    if args.flip_radial:
        radial_score = -radial_score

    longitudinal01 = normalize01(longitudinal_score)
    radial01 = normalize01(radial_score)

    longitudinal_img = args.digitize_pad + args.digitize_scale * longitudinal01
    radial_img = args.digitize_pad + args.digitize_scale * radial01

    # Spateo samples the generated field as field[row, col], while OpenCV contours
    # are represented as (x=col, y=row). Keep separate arrays for both conventions.
    spatial_digitize = np.column_stack([radial_img, longitudinal_img])
    contour_digitize = np.column_stack([longitudinal_img, radial_img])

    return {
        "longitudinal_axis": longitudinal_axis,
        "radial_axis": radial_axis,
        "longitudinal_score": longitudinal_score,
        "radial_score": radial_score,
        "longitudinal01": longitudinal01,
        "radial01": radial01,
        "spatial_digitize": spatial_digitize,
        "contour_digitize": contour_digitize,
    }


def build_convex_hull_contour(spatial_digitize, contour_digitize, digitize_pad):
    img_size = int(np.ceil(max(spatial_digitize.max(), contour_digitize.max()))) + digitize_pad + 1
    roi_img = np.zeros((img_size, img_size), dtype=np.uint8)

    hull_points = np.round(contour_digitize).astype(np.int32).reshape(-1, 1, 2)
    hull = cv2.convexHull(hull_points)
    cv2.drawContours(roi_img, [hull], contourIdx=0, color=255, thickness=cv2.FILLED)

    contours, _ = cv2.findContours(
        roi_img,
        cv2.RETR_EXTERNAL,
        cv2.CHAIN_APPROX_NONE,
    )

    if len(contours) == 0:
        raise RuntimeError("No ROI contour found.")

    areas = np.array([cv2.contourArea(c) for c in contours])
    ctr_idx = int(np.argmax(areas))

    # st.dd.digitize sizes its PDE field from max(adata.obsm[spatial_key]).
    # Keep contour vertices inside the field Spateo will allocate.
    field_rowmax = int(np.floor(spatial_digitize[:, 0].max()))
    field_colmax = int(np.floor(spatial_digitize[:, 1].max()))
    contours = tuple(
        np.stack(
            [
                np.clip(c[:, :, 0], 0, field_colmax),
                np.clip(c[:, :, 1], 0, field_rowmax),
            ],
            axis=2,
        ).astype(np.int32)
        for c in contours
    )

    return contours, ctr_idx, areas


def auto_corner_rows(longitudinal01, radial01):
    # Match the Spateo tutorial convention:
    # pnt_xy: min layer, min column
    # pnt_xY: max layer, min column
    # pnt_Xy: min layer, max column
    # pnt_XY: max layer, max column
    targets = {
        "pnt_xy": (0, 0),
        "pnt_xY": (0, 1),
        "pnt_Xy": (1, 0),
        "pnt_XY": (1, 1),
    }
    rows = {}
    for name, (column_target, layer_target) in targets.items():
        d2 = (longitudinal01 - column_target) ** 2 + (radial01 - layer_target) ** 2
        rows[name] = int(np.argmin(d2))
    return rows


def snap_corner_rows_to_contour(corner_rows, contour, contour_digitize, n_spots):
    for name, row_idx in corner_rows.items():
        if row_idx < 0 or row_idx >= n_spots:
            raise IndexError(f"{name} row index {row_idx} is outside xy file with {n_spots} rows")

    corner_points = {}
    corner_anchors = {}
    corner_snap_distances = {}

    for name, row_idx in corner_rows.items():
        anchor = contour_digitize[row_idx]
        d2_contour = (contour[:, 0] - anchor[0]) ** 2 + (contour[:, 1] - anchor[1]) ** 2
        contour_idx = int(np.argmin(d2_contour))
        corner_points[name] = tuple(map(int, contour[contour_idx]))
        corner_anchors[name] = anchor
        corner_snap_distances[name] = float(np.sqrt(d2_contour[contour_idx]))

    return corner_points, corner_anchors, corner_snap_distances


def plot_corner_qc(coords, contour_digitize, contour, corner_rows, corner_points, corner_anchors, plot_file):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    axes[0].scatter(contour_digitize[:, 0], contour_digitize[:, 1], s=8, alpha=0.5)
    axes[0].plot(contour[:, 0], contour[:, 1], color="black", linewidth=1)
    for name, row_idx in corner_rows.items():
        anchor = corner_anchors[name]
        point = corner_points[name]
        axes[0].scatter(anchor[0], anchor[1], s=70, marker="x", color="tab:red")
        axes[0].scatter(point[0], point[1], s=80)
        axes[0].plot([anchor[0], point[0]], [anchor[1], point[1]], color="tab:red", linewidth=0.8)
        axes[0].text(point[0], point[1], name, fontsize=9)
    axes[0].set_title("Granule corner rows snapped to ROI contour")
    axes[0].set_aspect("equal")

    axes[1].scatter(coords[:, 0], coords[:, 1], s=8, alpha=0.5)
    for name, row_idx in corner_rows.items():
        axes[1].scatter(coords[row_idx, 0], coords[row_idx, 1], s=80)
        axes[1].text(coords[row_idx, 0], coords[row_idx, 1], f"{name} ({row_idx})", fontsize=9)
    axes[1].set_title("Corner rows in original xy coordinates")
    axes[1].set_aspect("equal")

    fig.savefig(plot_file, dpi=180)
    plt.close(fig)


def plot_digitized_qc(coords, adata, plot_file):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    sc0 = axes[0].scatter(
        coords[:, 0],
        coords[:, 1],
        c=adata.obs["spateo_longitudinal"],
        s=10,
        cmap="viridis",
    )
    axes[0].set_title("Spateo longitudinal coordinate")
    axes[0].set_aspect("equal")
    fig.colorbar(sc0, ax=axes[0], fraction=0.046, pad=0.04)

    sc1 = axes[1].scatter(
        coords[:, 0],
        coords[:, 1],
        c=adata.obs["spateo_radial"],
        s=10,
        cmap="magma",
    )
    axes[1].set_title("Spateo radial coordinate")
    axes[1].set_aspect("equal")
    fig.colorbar(sc1, ax=axes[1], fraction=0.046, pad=0.04)

    fig.savefig(plot_file, dpi=180)
    plt.close(fig)


def main():
    args = parse_args()
    xy_file = args.xy_file or args.data_dir / "granule_xy.csv"
    output_file = args.output_file or args.data_dir / "granule_xy_digitized.csv"
    corner_rows_override = parse_corner_rows(args.corner_rows)

    args.plot_dir.mkdir(parents=True, exist_ok=True)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    xy, coords, adata = load_coordinates(xy_file)
    image_coords = pca_image_coordinates(coords, args)

    adata.obsm["spatial_digitize"] = image_coords["spatial_digitize"]
    adata.obs["pca_longitudinal_score"] = image_coords["longitudinal_score"]
    adata.obs["pca_radial_score"] = image_coords["radial_score"]

    contours, ctr_idx, areas = build_convex_hull_contour(
        image_coords["spatial_digitize"],
        image_coords["contour_digitize"],
        args.digitize_pad,
    )
    contour = contours[ctr_idx][:, 0, :]

    corner_rows = corner_rows_override or auto_corner_rows(
        image_coords["longitudinal01"],
        image_coords["radial01"],
    )
    corner_points, corner_anchors, corner_snap_distances = snap_corner_rows_to_contour(
        corner_rows,
        contour,
        image_coords["contour_digitize"],
        coords.shape[0],
    )

    print(f"Loaded {coords.shape[0]} granule spots from {xy_file}")
    print("longitudinal_axis:", image_coords["longitudinal_axis"])
    print("radial_axis:", image_coords["radial_axis"])
    print(
        "spatial_digitize row/col min/max:",
        image_coords["spatial_digitize"].min(axis=0),
        image_coords["spatial_digitize"].max(axis=0),
    )
    print(
        "contour_digitize x/y min/max:",
        image_coords["contour_digitize"].min(axis=0),
        image_coords["contour_digitize"].max(axis=0),
    )
    print(f"Found {len(contours)} contours; using ctr_idx={ctr_idx}, area={areas[ctr_idx]:.1f}")
    print("Corner rows and snapped Spateo contour points:")
    for name in CORNER_NAMES:
        idx = corner_rows[name]
        anchor = corner_anchors[name]
        point = corner_points[name]
        dist = corner_snap_distances[name]
        print(
            f"  {name}: granule_xy row={idx}, spot={adata.obs_names[idx]}, "
            f"anchor=({anchor[0]:.2f}, {anchor[1]:.2f}), "
            f"contour_point={point}, snap_distance={dist:.2f}"
        )

    if not args.no_plots:
        plot_corner_qc(
            coords,
            image_coords["contour_digitize"],
            contour,
            corner_rows,
            corner_points,
            corner_anchors,
            args.plot_dir / f"{args.plot_prefix}_corners.png",
        )

    st.dd.digitize(
        adata=adata,
        ctrs=contours,
        ctr_idx=ctr_idx,
        pnt_xy=corner_points["pnt_xy"],
        pnt_xY=corner_points["pnt_xY"],
        pnt_Xy=corner_points["pnt_Xy"],
        pnt_XY=corner_points["pnt_XY"],
        spatial_key="spatial_digitize",
        dgl_layer_key="spateo_radial_raw",
        dgl_column_key="spateo_longitudinal_raw",
        max_itr=args.max_itr,
        lh=args.lh,
        hh=args.hh,
    )

    adata.obs["spateo_radial"] = normalize01(adata.obs["spateo_radial_raw"].astype(float).to_numpy())
    adata.obs["spateo_longitudinal"] = normalize01(
        adata.obs["spateo_longitudinal_raw"].astype(float).to_numpy()
    )

    digitized_coords = pd.DataFrame(
        {
            "xcoord": coords[:, 0],
            "ycoord": coords[:, 1],
            "spateo_longitudinal": adata.obs["spateo_longitudinal"].to_numpy(),
            "spateo_radial": adata.obs["spateo_radial"].to_numpy(),
            "spateo_longitudinal_raw": adata.obs["spateo_longitudinal_raw"].to_numpy(),
            "spateo_radial_raw": adata.obs["spateo_radial_raw"].to_numpy(),
            "pca_longitudinal_score": image_coords["longitudinal_score"],
            "pca_radial_score": image_coords["radial_score"],
            "corner_role": "",
        },
        index=adata.obs_names,
    )
    for name, row_idx in corner_rows.items():
        digitized_coords.iloc[row_idx, digitized_coords.columns.get_loc("corner_role")] = name

    for key in ["spateo_longitudinal_raw", "spateo_radial_raw", "spateo_longitudinal", "spateo_radial"]:
        values = adata.obs[key].astype(float).to_numpy()
        print(f"{key}: zeros={(values == 0).sum()}, min={values.min():.4f}, max={values.max():.4f}")

    digitized_coords.to_csv(output_file, index=True)
    print(f"Wrote {output_file.resolve()}")

    if not args.no_plots:
        plot_digitized_qc(
            coords,
            adata,
            args.plot_dir / f"{args.plot_prefix}_coordinates.png",
        )
        print(f"Wrote QC plots to {args.plot_dir.resolve()}")


if __name__ == "__main__":
    main()
