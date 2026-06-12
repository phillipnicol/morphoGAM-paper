#!/usr/bin/env python

import argparse
import os
from pathlib import Path

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
        description="Run GASTON on the CA3 simulation grid."
    )
    parser.add_argument("--data-dir", type=Path, default=Path("../data"))
    parser.add_argument("--coords-file", type=Path, default=None)
    parser.add_argument("--output-file", type=Path, default=None)
    parser.add_argument("--gaston-output-root", type=Path, default=None)
    parser.add_argument("--patterns", nargs="+", default=list(DEFAULT_PATTERNS))
    parser.add_argument("--peak-t", nargs="+", type=float, default=list(DEFAULT_PEAK_T))
    parser.add_argument("--sd-r", nargs="+", type=float, default=list(DEFAULT_SD_R))
    parser.add_argument(
        "--replicates", nargs="+", type=int, default=list(DEFAULT_REPLICATES)
    )
    parser.add_argument("--num-dims", type=int, default=5)
    parser.add_argument("--glmpca-penalty", type=float, default=20)
    parser.add_argument("--glmpca-num-iters", type=int, default=30)
    parser.add_argument("--glmpca-eps", type=float, default=1e-4)
    parser.add_argument("--epochs", type=int, default=10000)
    parser.add_argument("--checkpoint", type=int, default=500)
    parser.add_argument("--num-restarts", type=int, default=2)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--isodepth-hidden", nargs="+", type=int, default=[20, 20])
    parser.add_argument("--expression-hidden", nargs="+", type=int, default=[20, 20])
    parser.add_argument("--optimizer", default="adam")
    parser.add_argument("--num-layers", type=int, default=1)
    parser.add_argument("--tutorial-threshold", type=float, default=0.1)
    parser.add_argument("--umi-threshold", type=float, default=0)
    parser.add_argument("--zero-fit-threshold", type=float, default=0)
    parser.add_argument("--isodepth-mult-factor", type=float, default=0.1)
    parser.add_argument(
        "--use-processed-inputs",
        action="store_true",
        help="Use A_processed/S_processed from process_files for isodepth labels. By default, mimic gaston_test.ipynb and use raw A/S.",
    )
    parser.add_argument(
        "--no-flip-isodepth",
        action="store_true",
        help="Do not flip isodepth after get_isodepth_labels. Rank by abs(slope) is unchanged by the flip.",
    )
    parser.add_argument(
        "--force-retrain",
        action="store_true",
        help="Train GASTON even when an output directory already exists.",
    )
    parser.add_argument("--max-datasets", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--rerun-failed",
        action="store_true",
        help="When resuming, rerun rows that have a nonempty gaston_error.",
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
        / f"ca3_simulation_{pattern}"
        / (
            "ca3_simulated_dataset_"
            f"peakt_{format_label(peak_t_value)}_"
            f"sd.r_{format_label(sd_r_value)}_"
            f"rep_{replicate:02d}.csv"
        )
    )


def gaston_output_dir(output_root, pattern, peak_t_value, sd_r_value, replicate):
    return output_root / (
        f"ca3_gaston_outputs_pattern_{pattern}_"
        f"peakt_{format_label(peak_t_value)}_"
        f"sd.r_{format_label(sd_r_value)}_"
        f"rep_{replicate:02d}"
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
    rank_ok = rows["gaston_rank"].notna() if "gaston_rank" in rows else False
    error_ok = (
        rows["gaston_error"].fillna("").astype(str).eq("")
        if "gaston_error" in rows
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


def empty_result(row, input_file, output_dir):
    return {
        "pattern": row["pattern"],
        "peak.t.value": row["peak.t.value"],
        "sd.r.value": row["sd.r.value"],
        "replicate": row["replicate"],
        "input_file": str(input_file),
        "gaston_output_dir": str(output_dir),
        "gaston_rank": np.nan,
        "gaston_slope": np.nan,
        "gaston_abs_slope": np.nan,
        "gaston_slope_tutorial_thresholded": np.nan,
        "gaston_p_value": np.nan,
        "gaston_slope_per_isodepth_sd": np.nan,
        "gaston_log_fc_across_isodepth_range": np.nan,
        "gaston_isodepth_sd": np.nan,
        "gaston_isodepth_range": np.nan,
        "gaston_error": "",
    }


def process_gaston_output(output_dir, process_nn_output):
    return process_nn_output.process_files(str(output_dir))


def train_gaston(S, A, output_dir, args, neural_net, torch):
    output_dir.mkdir(parents=True, exist_ok=True)
    S_torch, A_torch = neural_net.load_rescale_input_data(S, A)

    for restart in range(args.num_restarts):
        seed = args.seed + restart
        print(f"training GASTON restart {restart} with seed {seed}")
        out_dir_seed = output_dir / f"rep{restart}"
        out_dir_seed.mkdir(parents=True, exist_ok=True)

        torch.manual_seed(seed)
        neural_net.train(
            S_torch,
            A_torch,
            S_hidden_list=args.isodepth_hidden,
            A_hidden_list=args.expression_hidden,
            epochs=args.epochs,
            checkpoint=args.checkpoint,
            save_dir=str(out_dir_seed),
            optim=args.optimizer,
            seed=seed,
            save_final=True,
        )


def fit_glmpca(counts_mat, args):
    from glmpca import glmpca

    counts_mat_spot_by_gene = counts_mat.T
    num_genes = counts_mat.shape[0]
    top_gene_ix = np.argsort(np.sum(counts_mat_spot_by_gene, axis=0))[-num_genes:]
    counts_mat_glmpca = counts_mat_spot_by_gene[:, top_gene_ix]

    glmpca_res = glmpca.glmpca(
        counts_mat_glmpca.T,
        args.num_dims,
        fam="poi",
        penalty=args.glmpca_penalty,
        verbose=True,
        ctl={
            "maxIter": args.glmpca_num_iters,
            "eps": args.glmpca_eps,
            "optimizeTheta": True,
        },
    )
    return glmpca_res["factors"]


def run_gaston_for_dataset(counts_file, coords_file, output_dir, args):
    import torch
    from gaston import dp_related, neural_net, process_NN_output, segmented_fit

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    counts_df = pd.read_csv(counts_file, index_col=0)
    coords_df = pd.read_csv(coords_file)

    counts_mat = counts_df.to_numpy(dtype=np.int64)
    coords_mat = coords_df[["xcoord", "ycoord"]].to_numpy(dtype=np.float64)
    gene_labels = counts_df.index.astype(str).to_numpy()

    if counts_mat.shape[1] != coords_mat.shape[0]:
        raise ValueError(
            "Mismatched spots: "
            f"counts={counts_mat.shape[1]}, coords={coords_mat.shape[0]}"
        )

    A = fit_glmpca(counts_mat, args)
    S = coords_mat

    try_existing = output_dir.exists() and not args.force_retrain
    if try_existing:
        try:
            gaston_model, A_processed, S_processed = process_gaston_output(
                output_dir, process_NN_output
            )
        except Exception as exc:
            print(f"Could not reuse existing GASTON output {output_dir}: {exc}")
            train_gaston(S, A, output_dir, args, neural_net, torch)
            gaston_model, A_processed, S_processed = process_gaston_output(
                output_dir, process_NN_output
            )
    else:
        train_gaston(S, A, output_dir, args, neural_net, torch)
        gaston_model, A_processed, S_processed = process_gaston_output(
            output_dir, process_NN_output
        )

    A_for_labels = A_processed if args.use_processed_inputs else A
    S_for_labels = S_processed if args.use_processed_inputs else S
    gaston_isodepth, gaston_labels = dp_related.get_isodepth_labels(
        gaston_model,
        A_for_labels,
        S_for_labels,
        args.num_layers,
    )

    if not args.no_flip_isodepth:
        gaston_isodepth = np.max(gaston_isodepth) - gaston_isodepth
        gaston_labels = (args.num_layers - 1) - gaston_labels

    if len(np.unique(gaston_labels)) != args.num_layers:
        raise ValueError(
            f"Expected {args.num_layers} GASTON layer, found {np.unique(gaston_labels)}."
        )

    counts_mat_for_pw = counts_mat
    if counts_mat_for_pw.shape[0] != len(gaston_isodepth):
        counts_mat_for_pw = counts_mat_for_pw.T

    if counts_mat_for_pw.shape[0] != len(gaston_isodepth):
        raise ValueError(
            "Counts rows must match isodepth length: "
            f"counts={counts_mat_for_pw.shape[0]}, isodepth={len(gaston_isodepth)}"
        )

    idx_kept = np.arange(counts_mat_for_pw.shape[1])
    pw_fit_dict = segmented_fit.pw_linear_fit(
        counts_mat_for_pw,
        gaston_labels,
        gaston_isodepth,
        None,
        [],
        idx_kept=idx_kept,
        zero_fit_threshold=args.zero_fit_threshold,
        t=np.inf,
        umi_threshold=args.umi_threshold,
        isodepth_mult_factor=args.isodepth_mult_factor,
    )

    slope_mat, intercept_mat, discont_mat, pv_mat = pw_fit_dict["all_cell_types"]
    if slope_mat.shape[1] != 1:
        raise ValueError(
            f"Expected one slope column for one layer, found {slope_mat.shape[1]}."
        )

    slope = slope_mat[:, 0]
    p_values = pv_mat[:, 0]
    slope_thresholded = np.where(p_values < args.tutorial_threshold, slope, 0.0)
    isodepth_sd = np.std(gaston_isodepth)
    isodepth_range = np.ptp(gaston_isodepth)

    gene_slopes = pd.DataFrame(
        {
            "gene": gene_labels[idx_kept],
            "slope": slope,
            "slope_tutorial_thresholded": slope_thresholded,
            "slope_per_isodepth_sd": slope * isodepth_sd,
            "log_fc_across_isodepth_range": slope * isodepth_range,
            "p_value": p_values,
        }
    )

    gene_slopes_ranked = (
        gene_slopes.assign(abs_slope=gene_slopes["slope"].abs())
        .sort_values("abs_slope", ascending=False, kind="mergesort")
        .reset_index(drop=True)
    )
    gene_slopes_ranked["abs_slope_rank"] = np.arange(
        1, len(gene_slopes_ranked) + 1
    )

    gene_1_rows = gene_slopes_ranked[gene_slopes_ranked["gene"].astype(str) == "1"]
    if gene_1_rows.empty:
        raise ValueError("Gene '1' was not found in GASTON slope results.")

    gene_1 = gene_1_rows.iloc[0]
    return {
        "gaston_rank": int(gene_1["abs_slope_rank"]),
        "gaston_slope": float(gene_1["slope"]),
        "gaston_abs_slope": float(gene_1["abs_slope"]),
        "gaston_slope_tutorial_thresholded": float(
            gene_1["slope_tutorial_thresholded"]
        ),
        "gaston_p_value": float(gene_1["p_value"]),
        "gaston_slope_per_isodepth_sd": float(gene_1["slope_per_isodepth_sd"]),
        "gaston_log_fc_across_isodepth_range": float(
            gene_1["log_fc_across_isodepth_range"]
        ),
        "gaston_isodepth_sd": float(isodepth_sd),
        "gaston_isodepth_range": float(isodepth_range),
    }


def main():
    args = parse_args()
    data_dir = args.data_dir
    coords_file = args.coords_file or data_dir / "ca3_xy.csv"
    output_file = args.output_file or data_dir / "ca3_gaston_model_results.csv"
    output_root = args.gaston_output_root or data_dir / "gaston_outputs"

    grid = simulation_grid(args)
    if args.dry_run:
        for row in grid:
            input_file = dataset_path(
                data_dir,
                row["pattern"],
                row["peak.t.value"],
                row["sd.r.value"],
                row["replicate"],
            )
            out_dir = gaston_output_dir(
                output_root,
                row["pattern"],
                row["peak.t.value"],
                row["sd.r.value"],
                row["replicate"],
            )
            print(f"{input_file} -> {out_dir}")
        print(f"Would process {len(grid)} datasets.")
        return

    output_file.parent.mkdir(parents=True, exist_ok=True)
    if args.overwrite and output_file.exists():
        output_file.unlink()
        print(f"Removed existing GASTON results file because --overwrite was set: {output_file}")

    if output_file.exists() and not args.overwrite:
        all_results = pd.read_csv(output_file)
        all_results = keep_current_grid(all_results, grid)
    else:
        all_results = pd.DataFrame()

    done = existing_keys(all_results, args.rerun_failed)
    output_root.mkdir(parents=True, exist_ok=True)

    for i, row in enumerate(grid, start=1):
        key = key_tuple(row)
        if key in done:
            print(f"[{i}/{len(grid)}] Skipping existing GASTON result: {key}")
            continue

        input_file = dataset_path(
            data_dir,
            row["pattern"],
            row["peak.t.value"],
            row["sd.r.value"],
            row["replicate"],
        )
        out_dir = gaston_output_dir(
            output_root,
            row["pattern"],
            row["peak.t.value"],
            row["sd.r.value"],
            row["replicate"],
        )
        print(f"[{i}/{len(grid)}] Running GASTON: {input_file}")
        result = empty_result(row, input_file, out_dir)

        try:
            result.update(run_gaston_for_dataset(input_file, coords_file, out_dir, args))
        except Exception as exc:
            result["gaston_error"] = str(exc)

        all_results = without_key(all_results, key)
        all_results = pd.concat([all_results, pd.DataFrame([result])], ignore_index=True)
        all_results.to_csv(output_file, index=False)
        print(f"Checkpointed GASTON results to {output_file}")

    print(f"Wrote GASTON results to {output_file.resolve()}")


if __name__ == "__main__":
    main()
