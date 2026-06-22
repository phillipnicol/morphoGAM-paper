#!/usr/bin/env python

import argparse
import os
import warnings
from pathlib import Path

os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/cache-spateo")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mpl-spateo")
os.environ.setdefault("NUMBA_CACHE_DIR", "/private/tmp/numba-spateo")

import numpy as np
import pandas as pd


DEFAULT_PATTERNS = ("peak", "gradient", "nonparam")
DEFAULT_PEAK_T = (1.0,)
DEFAULT_SD_R = (0.5, 1.0, 2.0)
DEFAULT_REPLICATES = tuple(range(1, 11))
KEY_COLUMNS = ("pattern", "peak.t.value", "sd.r.value", "replicate")


def format_label(x):
    return f"{x:g}".replace(".", "p")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Spateo GLM on the granule simulation grid."
    )
    parser.add_argument("--data-dir", type=Path, default=Path("../data"))
    parser.add_argument(
        "--digitized-coords",
        type=Path,
        default=None,
        help="CSV containing spateo_longitudinal and spateo_radial. Defaults to DATA_DIR/granule_xy_digitized.csv.",
    )
    parser.add_argument("--output-file", type=Path, default=None)
    parser.add_argument("--patterns", nargs="+", default=list(DEFAULT_PATTERNS))
    parser.add_argument("--peak-t", nargs="+", type=float, default=list(DEFAULT_PEAK_T))
    parser.add_argument("--sd-r", nargs="+", type=float, default=list(DEFAULT_SD_R))
    parser.add_argument(
        "--replicates", nargs="+", type=int, default=list(DEFAULT_REPLICATES)
    )
    parser.add_argument(
        "--full-model", default="~cr(t, df=9) + cr(r, df=5)"
    )
    parser.add_argument("--reduced-model", default="~cr(r, df=5)")
    parser.add_argument("--max-datasets", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--rerun-failed",
        action="store_true",
        help="When resuming, rerun rows that have a nonempty spateo_error.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print the datasets that would be processed.",
    )
    return parser.parse_args()


def simulation_grid(args):
    rows = []
    for pattern in args.patterns:
        for peak_t_value in args.peak_t:
            for sd_r_value in args.sd_r:
                for replicate in args.replicates:
                    rows.append(
                        {
                            "pattern": pattern,
                            "peak.t.value": peak_t_value,
                            "sd.r.value": sd_r_value,
                            "replicate": replicate,
                        }
                    )

    if args.max_datasets is not None:
        rows = rows[: args.max_datasets]
    return rows


def dataset_path(data_dir, pattern, peak_t_value, sd_r_value, replicate):
    return (
        data_dir
        / f"granule_simulation_{pattern}"
        / (
            "granule_simulated_dataset_"
            f"peakt_{format_label(peak_t_value)}_"
            f"sd.r_{format_label(sd_r_value)}_"
            f"rep_{replicate:02d}.csv"
        )
    )


def key_tuple(row):
    return (
        str(row["pattern"]),
        float(row["peak.t.value"]),
        float(row["sd.r.value"]),
        int(row["replicate"]),
    )


def existing_keys(existing_results, rerun_failed):
    if existing_results.empty:
        return set()

    rows = existing_results
    rank_ok = rows["spateo_rank"].notna() if "spateo_rank" in rows else False
    error_ok = (
        rows["spateo_error"].fillna("").astype(str).eq("")
        if "spateo_error" in rows
        else True
    )
    rows = rows[rank_ok & error_ok]

    return {key_tuple(row) for _, row in rows.iterrows()}


def keep_current_grid(results, grid):
    if results.empty:
        return results

    grid_keys = {key_tuple(row) for row in grid}
    keep = [key_tuple(row) in grid_keys for _, row in results.iterrows()]
    return results.loc[keep].copy()


def without_key(results, key):
    if results.empty:
        return results

    mask = (
        (results["pattern"].astype(str) == key[0])
        & (results["peak.t.value"].astype(float) == key[1])
        & (results["sd.r.value"].astype(float) == key[2])
        & (results["replicate"].astype(int) == key[3])
    )
    return results.loc[~mask].copy()


def empty_result(row, input_file):
    return {
        "pattern": row["pattern"],
        "peak.t.value": row["peak.t.value"],
        "sd.r.value": row["sd.r.value"],
        "replicate": row["replicate"],
        "input_file": str(input_file),
        "spateo_rank": np.nan,
        "spateo_qval": np.nan,
        "spateo_pval": np.nan,
        "spateo_log_likelihood": np.nan,
        "spateo_status": "",
        "spateo_family": "",
        "spateo_error": "",
    }


def numeric_value(row, column):
    if column not in row.index:
        return np.nan
    value = row[column]
    return float(value) if pd.notna(value) else np.nan


def run_spateo_for_dataset(counts_file, digitized_path, full_model, reduced_model):
    import anndata as ad
    from spateo.tools.glm import glm_degs

    counts_df = pd.read_csv(counts_file, index_col=0)
    digitized_coords = pd.read_csv(digitized_path, index_col=0)

    required_cols = ["xcoord", "ycoord", "spateo_longitudinal", "spateo_radial"]
    missing_cols = [col for col in required_cols if col not in digitized_coords.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in {digitized_path}: {missing_cols}")

    raw_counts = counts_df.to_numpy(dtype=float).T
    gene_names = counts_df.index.astype(str).to_list()

    if raw_counts.shape[0] != digitized_coords.shape[0]:
        raise ValueError(
            "Mismatched spots: "
            f"counts={raw_counts.shape[0]}, digitized_coords={digitized_coords.shape[0]}"
        )

    library_size = raw_counts.sum(axis=1)
    positive_library = library_size[library_size > 0]
    median_library_size = np.median(positive_library) if len(positive_library) else 1.0
    size_factors = np.where(library_size > 0, library_size / median_library_size, 1.0)
    normalized_counts = raw_counts / size_factors[:, None]

    obs = pd.DataFrame(
        {
            "xcoord": digitized_coords["xcoord"].to_numpy(),
            "ycoord": digitized_coords["ycoord"].to_numpy(),
            "t": digitized_coords["spateo_longitudinal"].to_numpy(),
            "r": digitized_coords["spateo_radial"].to_numpy(),
            "size_factor": size_factors,
            "library_size": library_size,
        },
        index=[f"spot_{i + 1}" for i in range(raw_counts.shape[0])],
    )

    adata = ad.AnnData(
        X=normalized_counts,
        obs=obs,
        var=pd.DataFrame(index=gene_names),
    )
    adata.layers["counts"] = raw_counts
    adata.obsm["spatial"] = digitized_coords[["xcoord", "ycoord"]].to_numpy()
    adata.var_names_make_unique()

    glm_key = "glm_degs_t_linear"
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Negative binomial dispersion parameter alpha not set.*",
        )
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        glm_degs(
            adata,
            key_added=glm_key,
            fullModelFormulaStr=full_model,
            reducedModelFormulaStr=reduced_model,
            qval_threshold=None,
            llf_threshold=None,
            inplace=True,
        )

    glm_hits = adata.uns[glm_key]["glm_result"].copy()
    glm_hits.index = glm_hits.index.astype(str)
    sort_cols = [col for col in ["qval", "pval", "log-likelihood"] if col in glm_hits]
    glm_hits = glm_hits.sort_values(sort_cols, ascending=True, na_position="last")

    gene_positions = np.flatnonzero(glm_hits.index.to_numpy(dtype=str) == "1")
    if len(gene_positions) == 0:
        raise ValueError("Gene '1' was not found in Spateo GLM results.")

    gene_1_rank = int(gene_positions[0] + 1)
    gene_1 = glm_hits.iloc[gene_positions[0]]

    return {
        "spateo_rank": gene_1_rank,
        "spateo_qval": numeric_value(gene_1, "qval"),
        "spateo_pval": numeric_value(gene_1, "pval"),
        "spateo_log_likelihood": numeric_value(gene_1, "log-likelihood"),
        "spateo_status": str(gene_1.get("status", "")),
        "spateo_family": str(gene_1.get("family", "")),
    }


def main():
    args = parse_args()
    data_dir = args.data_dir
    digitized_path = args.digitized_coords or data_dir / "granule_xy_digitized.csv"
    output_file = args.output_file or data_dir / "granule_spateo_model_results.csv"

    grid = simulation_grid(args)
    if args.dry_run:
        for row in grid:
            print(
                dataset_path(
                    data_dir,
                    row["pattern"],
                    row["peak.t.value"],
                    row["sd.r.value"],
                    row["replicate"],
                )
            )
        print(f"Would process {len(grid)} datasets.")
        return

    output_file.parent.mkdir(parents=True, exist_ok=True)
    if args.overwrite and output_file.exists():
        output_file.unlink()
        print(f"Removed existing Spateo results file because --overwrite was set: {output_file}")

    if output_file.exists() and not args.overwrite:
        all_results = pd.read_csv(output_file)
        all_results = keep_current_grid(all_results, grid)
    else:
        all_results = pd.DataFrame()

    done = existing_keys(all_results, args.rerun_failed)

    for i, row in enumerate(grid, start=1):
        key = key_tuple(row)
        if key in done:
            print(f"[{i}/{len(grid)}] Skipping existing Spateo result: {key}")
            continue

        input_file = dataset_path(
            data_dir,
            row["pattern"],
            row["peak.t.value"],
            row["sd.r.value"],
            row["replicate"],
        )
        print(f"[{i}/{len(grid)}] Running Spateo: {input_file}")
        result = empty_result(row, input_file)

        try:
            result.update(
                run_spateo_for_dataset(
                    input_file,
                    digitized_path,
                    args.full_model,
                    args.reduced_model,
                )
            )
        except Exception as exc:
            result["spateo_error"] = str(exc)

        all_results = without_key(all_results, key)
        all_results = pd.concat([all_results, pd.DataFrame([result])], ignore_index=True)
        all_results.to_csv(output_file, index=False)
        print(f"Checkpointed Spateo results to {output_file}")

    print(f"Wrote Spateo results to {output_file.resolve()}")


if __name__ == "__main__":
    main()
