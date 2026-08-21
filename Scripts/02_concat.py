#!/path/to/envs/scanpy/bin/python
import sys, argparse, warnings, platform, os, glob, time
from datetime import datetime
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


def plot_qc_metrics(qc_rows, outpath, prefix, dpi=300, stage="before"):
    """Plot QC violin + scatter from collected per-slice QC metrics.

    Writes {prefix}{stage}_QC_violin.png/pdf, {prefix}{stage}_QC_scatter.png/pdf
    (stage = 'before' | 'after') — filenames must match Report_config.
    """
    df = pd.concat(qc_rows, ignore_index=True)
    df['slice_id'] = df['slice_id'].astype('category')

    features = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']

    # --- violin plots (1x3 grid, seaborn) ---
    # cut=0: 截断 KDE 在数据范围边界处，避免平滑溢出到 0 以下
    # set_ylim(bottom=0): matplotlib 自动 ylim 会加 margin 使下界 < 0，
    #   n_genes / total_counts / pct_counts_mt 均不可能为负，强制下界=0
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
    for ax, feature in zip(axes.flatten(), features):
        sns.violinplot(data=df, x='slice_id', y=feature, ax=ax,
                       saturation=0.8, inner=None, cut=0)
        ax.set_ylim(bottom=0)
        ax.set_xlabel('slice_id')
        ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.savefig("{}/{}{}_QC_violin.png".format(outpath, prefix, stage), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}{}_QC_violin.pdf".format(outpath, prefix, stage), dpi=dpi, bbox_inches='tight')
    plt.close()

    # --- scatter plots (1x2 grid) ---
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14, 5))
    sns.scatterplot(data=df, x='total_counts', y='n_genes_by_counts', hue='slice_id',
                    s=5, alpha=0.6, ax=ax[0], legend=False)
    sns.scatterplot(data=df, x='total_counts', y='pct_counts_mt', hue='slice_id',
                    s=5, alpha=0.6, ax=ax[1], legend=False)
    plt.tight_layout()
    plt.savefig("{}/{}{}_QC_scatter.png".format(outpath, prefix, stage), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}{}_QC_scatter.pdf".format(outpath, prefix, stage), dpi=dpi, bbox_inches='tight')
    plt.close()


def concat_h5ad(input_dir, outpath, prefix,
                join_method, merge_method, uns_merge_method, label_name,
                q_low, q_high, min_cells, max_mt):
    """Per-sample QC + normalize + concatenate multiple H5AD files.

    For each slice:
      1. Per-slice QC: keep cells within quantile range of n_genes/total_counts
         and below the mt% threshold.
      2. Per-slice normalize: normalize_total(target_sum=1e4) + log1p.
      3. Collect into a list, then concatenate.
    """

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Find all H5AD files                                                  #
    # ------------------------------------------------------------------ #
    h5ad_files = sorted(glob.glob(os.path.join(input_dir, "*.h5ad")))

    if not h5ad_files:
        print(f"Error: No H5AD files found in {input_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(h5ad_files)} H5AD files:")
    for f in h5ad_files:
        print(f"  - {os.path.basename(f)}")

    # ------------------------------------------------------------------ #
    # Per-sample QC + normalize (loop over each slice)                    #
    # ------------------------------------------------------------------ #
    print(f"\n>>> Per-sample QC + normalize (quantiles [{q_low:.0%}, {q_high:.0%}], mt<={max_mt})...\n")
    adata_list = []
    sample_names = []
    threshold_table = []
    qc_before_rows = []   # collect QC metrics before filtering (for plots)
    qc_after_rows = []    # collect QC metrics after filtering (for plots)

    for i, h5ad_file in enumerate(h5ad_files):
        try:
            adata = ad.read_h5ad(h5ad_file)
        except Exception as e:
            print(f"  [ERROR] Failed to load {h5ad_file}: {str(e)}", file=sys.stderr)
            sys.exit(1)

        # Get sample name from slice_id if available, otherwise from filename
        if "slice_id" in adata.obs.columns:
            sample_name = adata.obs["slice_id"].iloc[0]
        elif "orig.ident" in adata.obs.columns:
            sample_name = adata.obs["orig.ident"].iloc[0]
        else:
            sample_name = os.path.splitext(os.path.basename(h5ad_file))[0]
        sample_names.append(sample_name)

        # Calculate QC metrics so n_genes_by_counts / total_counts / pct_counts_mt exist
        # 线粒体基因匹配：MT / mt 开头（含 MT- / mt- 横线形式），兼容不同物种命名习惯
        # （人类 MT-ND1，小鼠 mt-Nd1，鱼类 mt-nd1 等；LOC 前缀数据无 mt 基因时全为 False）
        adata.var['mt'] = adata.var_names.str.match(r'^(MT|mt)')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                                   log1p=False, inplace=True)

        # Record QC metrics BEFORE filtering (for before_QC plots)
        qc_before_rows.append(pd.DataFrame({
            'slice_id': sample_name,
            'n_genes_by_counts': adata.obs['n_genes_by_counts'].values,
            'total_counts': adata.obs['total_counts'].values,
            'pct_counts_mt': adata.obs['pct_counts_mt'].values,
        }))

        n_in = adata.n_obs
        n_genes_in = adata.n_vars

        # Per-sample QC (computed on raw counts, per slice independently)
        if n_in < 10:
            # Too few cells to estimate quantiles — keep all
            threshold_table.append((sample_name, n_in, n_in, '-', '-', '-', '-'))
            print(f"  [{i+1}/{len(h5ad_files)}] {sample_name}: n={n_in} (too few, kept all)")
        else:
            ng_q = np.quantile(adata.obs['n_genes_by_counts'], [q_low, q_high])
            tc_q = np.quantile(adata.obs['total_counts'],    [q_low, q_high])
            ng_lo, ng_hi = int(ng_q[0]), int(ng_q[1])
            tc_lo, tc_hi = int(tc_q[0]), int(tc_q[1])
            keep_mask = (
                (adata.obs['n_genes_by_counts'] >= ng_lo) &
                (adata.obs['n_genes_by_counts'] <= ng_hi) &
                (adata.obs['total_counts']      >= tc_lo) &
                (adata.obs['total_counts']      <= tc_hi) &
                (adata.obs['pct_counts_mt']     <= max_mt)
            ).values
            n_keep = int(keep_mask.sum())
            threshold_table.append((sample_name, n_in, n_keep, ng_lo, ng_hi, tc_lo, tc_hi))
            print(f"  [{i+1}/{len(h5ad_files)}] {sample_name}: n={n_in} → {n_keep} "
                  f"(ng=[{ng_lo},{ng_hi}], tc=[{tc_lo},{tc_hi}], mt<={max_mt})")

        # Apply filter + per-sample normalize
        sc.pp.filter_genes(adata, min_cells=min_cells)
        adata = adata[keep_mask].copy()

        # Record QC metrics AFTER filtering (for after_QC plots)
        qc_after_rows.append(pd.DataFrame({
            'slice_id': sample_name,
            'n_genes_by_counts': adata.obs['n_genes_by_counts'].values,
            'total_counts': adata.obs['total_counts'].values,
            'pct_counts_mt': adata.obs['pct_counts_mt'].values,
        }))

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        # Preserve log1p layer — required by downstream scripts (06_align reads
        # adata.layers['log1p'] to set X; V2 stored this layer in 03, V3 moved
        # normalization here so the layer must be created here).
        adata.layers['log1p'] = adata.X.copy()
        adata_list.append(adata)

    # ------------------------------------------------------------------ #
    # Print per-slice threshold table                                      #
    # ------------------------------------------------------------------ #
    print("\nPer-slice thresholds:")
    print("  {:<20s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}".format(
        "slice", "n_in", "n_keep", "ng_lo", "ng_hi", "tc_lo", "tc_hi"))
    for row in threshold_table:
        print("  {:<20s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}".format(*[str(x) for x in row]))
    print()

    # ------------------------------------------------------------------ #
    # Plot QC metrics before/after (filenames must match Report_config)   #
    # ------------------------------------------------------------------ #
    if qc_before_rows:
        plot_qc_metrics(qc_before_rows, outpath, prefix, dpi=300, stage="before")
    if qc_after_rows:
        plot_qc_metrics(qc_after_rows, outpath, prefix, dpi=300, stage="after")

    # ------------------------------------------------------------------ #
    # Concatenate H5AD files (X is now log1p-normalized per slice)        #
    # ------------------------------------------------------------------ #
    print(f">>> Concatenating {len(adata_list)} files...")
    print(f"    join      : {join_method}")
    print(f"    merge     : {merge_method}")
    print(f"    uns_merge : {uns_merge_method}")
    print(f"    label     : {label_name}")
    print(f"    keys      : {sample_names}")

    adata_concat = ad.concat(
        adata_list,
        join=join_method,
        merge=merge_method,
        uns_merge=uns_merge_method,
        label=label_name,
        keys=sample_names,
        index_unique="-"
    )

    print(f"\nConcatenated result:")
    print(f"    Cells: {adata_concat.n_obs}")
    print(f"    Genes: {adata_concat.n_vars}")
    print(f"    Batches: {adata_concat.obs[label_name].nunique()}")
    print(f"    X dtype : {adata_concat.X.dtype} (should be float = log1p-normalized)")

    # ------------------------------------------------------------------ #
    # Save concatenated H5AD file                                          #
    # ------------------------------------------------------------------ #
    output_name = f"{prefix}concated.h5ad"
    output_path = os.path.join(outpath, output_name)

    print(f"\n>>> Saving to: {output_path}")
    adata_concat.write_h5ad(output_path, compression="gzip")

    elapsed = time.time() - pipeline_start
    print(f"\nFile: {output_name} saved successfully !")
    print("Total time   : {:.1f} s\n".format(elapsed))

    return adata_concat


def main():
    parser = argparse.ArgumentParser(
        description="Per-sample QC + normalize + concatenate H5AD files into a single log1p-normalized matrix."
    )

    # Input/Output arguments
    parser.add_argument("-I", "--input",
                        type=str,
                        required=True,
                        metavar="INPUT_DIR",
                        help="Directory containing H5AD files to process")
    parser.add_argument("-O", "--outpath",
                        type=str,
                        default=".",
                        metavar="OUTPUT_DIR",
                        help="Output directory for the concatenated H5AD file (default: current directory)")
    parser.add_argument("-P", "--prefix",
                        type=str,
                        default="",
                        metavar="PREFIX",
                        help="Prefix for the output H5AD filename (default: none)")

    # Concatenation arguments
    parser.add_argument("--join",
                        type=str,
                        default="inner",
                        choices=["outer", "inner"],
                        metavar="JOIN",
                        help="Join method: 'outer' keeps all genes, 'inner' keeps only common genes "
                             "(default: inner). Inner avoids sparse-zero dilution when merging slices "
                             "with non-overlapping gene panels.")
    parser.add_argument("--merge",
                        type=str,
                        default="same",
                        choices=["same", "unique", "first", "outer", "none"],
                        metavar="MERGE",
                        help="Merge method for obs/var: 'same' keeps common columns, 'unique' keeps all unique columns (default: same)")
    parser.add_argument("--uns-merge",
                        type=str,
                        default="unique",
                        choices=["same", "unique", "first", "outer", "none"],
                        metavar="UNS_MERGE",
                        help="Merge method for uns dict keys: 'unique' keeps all unique keys (default: unique)")
    parser.add_argument("--label",
                        type=str,
                        default="batch",
                        metavar="LABEL",
                        help="Column name for batch labels in obs (default: batch)")

    # Per-sample QC arguments
    parser.add_argument("-qL", "--quantile_low",
                        type=float,
                        default=0.01,
                        metavar="FLOAT",
                        help="Per-slice lower quantile for n_genes / total_counts (default: 0.01)")
    parser.add_argument("-qH", "--quantile_high",
                        type=float,
                        default=0.99,
                        metavar="FLOAT",
                        help="Per-slice upper quantile for n_genes / total_counts (default: 0.99)")
    parser.add_argument("-minC", "--min_cells",
                        type=int,
                        default=3,
                        metavar="INT",
                        help="Minimum cells per gene (default: 3)")
    parser.add_argument("-maxMT", "--max_mt",
                        type=float,
                        default=10.0,
                        metavar="FLOAT",
                        help="Maximum mitochondrial percentage (default: 10.0)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    input_dir = args.input
    outpath = args.outpath
    prefix = args.prefix
    join_method = args.join
    merge_method = None if args.merge == "none" else args.merge
    uns_merge_method = None if args.uns_merge == "none" else args.uns_merge
    label_name = args.label
    q_low = args.quantile_low
    q_high = args.quantile_high
    min_cells = args.min_cells
    max_mt = args.max_mt

    # ------------------------------------------------------------------ #
    # Print run information                                                #
    # ------------------------------------------------------------------ #
    print("\n=========================== H5AD Per-sample QC + Concat ===========================\n")
    print("Current time   :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Operating system:", platform.system(), platform.release())
    print("Platform       :", platform.platform())
    print("Working directory:", os.getcwd())
    print("Input directory:", input_dir)
    print("Output directory:", outpath)
    print("Output prefix  :", f"'{prefix}'" if prefix else "(none)")
    print("Join method    :", join_method)
    print("Merge method   :", merge_method)
    print("Uns merge method:", uns_merge_method)
    print("Batch label    :", label_name)
    print(f"QC quantiles   : [{q_low:.0%}, {q_high:.0%}]")
    print(f"min_cells      : {min_cells}")
    print(f"max_mt         : {max_mt}")
    print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
    print("\n=================================================================================\n")

    # ------------------------------------------------------------------ #
    # Validate input directory                                             #
    # ------------------------------------------------------------------ #
    if not os.path.isdir(input_dir):
        print(f"Error: Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # Create output directory if needed                                    #
    # ------------------------------------------------------------------ #
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # ------------------------------------------------------------------ #
    # Per-sample QC + normalize + concatenate                             #
    # ------------------------------------------------------------------ #
    concat_h5ad(input_dir, outpath, prefix, join_method, merge_method, uns_merge_method, label_name,
                q_low, q_high, min_cells, max_mt)


if __name__ == "__main__":
    main()