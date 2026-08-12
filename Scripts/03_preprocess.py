#!/path/to/envs/scanpy/bin/python
import sys, argparse, warnings, os, logging, time
from datetime import datetime
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('fontTools').setLevel(logging.WARNING)


def main():
    parser = argparse.ArgumentParser(
        description="Single-cell RNA-seq analysis pipeline: per-sample QC -> "
                    "Normalization -> Feature selection -> PCA -> Batch correction -> "
                    "Clustering -> Marker gene identification"
    )

    # Input/Output arguments
    parser.add_argument("-I", "--input", required=True, metavar="H5AD_FILE",
                        help="Path to input h5ad file")
    parser.add_argument("-O", "--outpath", type=str, default=".", metavar="OUTPUT_DIR",
                        help="Output directory (default: current directory)")
    parser.add_argument("-P", "--prefix", type=str, default="", metavar="PREFIX",
                        help="Prefix for output files (default: none)")

    # QC arguments.
    #   n_genes / total_counts thresholds are derived per slice from quantiles,
    #   so no manual min/max is needed for those. Only global filters remain:
    #     -minC  : drop genes detected in fewer than N cells (per-slice aware via batch_key)
    #     -maxMT : drop cells with mitochondrial % above threshold
    parser.add_argument("-minC", "--min_cells", type=int, default=3, metavar="INT",
                        help="Minimum cells per gene (default: 3)")
    parser.add_argument("-maxMT", "--max_mt", type=float, default=10.0, metavar="FLOAT",
                        help="Maximum mitochondrial percentage (default: 10.0)")
    parser.add_argument("-qL", "--quantile_low", type=float, default=0.01, metavar="FLOAT",
                        help="Per-slice lower quantile for n_genes / total_counts (default: 0.01)")
    parser.add_argument("-qH", "--quantile_high", type=float, default=0.99, metavar="FLOAT",
                        help="Per-slice upper quantile for n_genes / total_counts (default: 0.99)")

    # Batch correction arguments
    parser.add_argument("-BK", "--batch_key", type=str, default="slice_id", metavar="STR",
                        help="Batch key in adata.obs for per-slice QC and batch correction "
                             "(default: slice_id)")

    # Clustering arguments
    parser.add_argument("-nPC", "--n_pcs", type=int, default=30, metavar="INT",
                        help="Number of PCs for clustering (default: 30)")
    parser.add_argument("-nHVG", "--n_top_genes", type=int, default=2000, metavar="INT",
                        help="Number of highly variable genes (default: 2000)")
    parser.add_argument("-res", "--resolution", type=float, default=1.0, metavar="FLOAT",
                        help="Leiden clustering resolution (default: 1.0)")

    # Plotting arguments
    parser.add_argument("-dpi", "--dpi", type=int, default=300, metavar="INT",
                        help="DPI for saved figures (default: 300)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    input_file    = args.input
    outpath       = args.outpath
    prefix        = args.prefix
    min_cells     = args.min_cells
    max_mt        = args.max_mt
    quantile_low  = args.quantile_low
    quantile_high = args.quantile_high
    batch_key     = args.batch_key
    n_pcs         = args.n_pcs
    n_top_genes   = args.n_top_genes
    resolution    = args.resolution
    dpi           = args.dpi

    # ------------------------------------------------------------------ #
    # Setup environment                                                    #
    # ------------------------------------------------------------------ #
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams["font.sans-serif"] = "Arial"
    sc.set_figure_params(dpi=80, dpi_save=dpi, figsize=(8, 8), facecolor="white")
    sc.settings.verbosity = 1
    os.makedirs(outpath, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Print run information                                                #
    # ------------------------------------------------------------------ #
    print("\n==================== scRNA-seq Analysis Pipeline ====================\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Working directory:", os.getcwd())
    print("Input file      :", input_file)
    print("Output directory:", outpath)
    print("Output prefix   :", prefix)
    print("Batch key       :", batch_key)
    print("Clustering      : n_pcs={}, n_top_genes={}, resolution={}".format(
          n_pcs, n_top_genes, resolution))
    print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
    print("\n======================================================================\n")

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Step 1: Load data                                                    #
    # ------------------------------------------------------------------ #
    print(">>> Step 1: Loading data\n")
    adata = sc.read_h5ad(input_file)
    print("Data shape: {} cells x {} genes\n".format(adata.n_obs, adata.n_vars))

    # ------------------------------------------------------------------ #
    # Step 2: Quality control                                              #
    # ------------------------------------------------------------------ #
    print(">>> Step 2: Quality control\n")

    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'],
                               percent_top=None, log1p=False, inplace=True)

    # Plot QC metrics before filtering
    features = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']

    if batch_key in adata.obs.columns:
        # Plot with batch grouping
        print("Plotting QC metrics grouped by: {}\n".format(batch_key))

        # Violin plots grouped by batch (1x3 grid, no stripplot)
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
        for ax, feature in zip(axes.flatten(), features):
            sc.pl.violin(adata, keys=feature, groupby=batch_key, stripplot=False,
                         rotation=45, ax=ax, show=False)
        plt.tight_layout()
        plt.savefig("{}/{}before_QC_violin.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}before_QC_violin.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

        # Scatter plots colored by batch (shared legend at bottom)
        n_batches = len(adata.obs[batch_key].astype('category').cat.categories)
        # Use enough colors for all batches
        if n_batches <= 20:
            colors = plt.cm.tab20(range(n_batches))
        elif n_batches <= 40:
            colors = np.vstack([plt.cm.tab20(range(20)), plt.cm.tab20b(range(n_batches - 20))])
        else:
            colors = plt.cm.hsv(np.linspace(0, 1, n_batches, endpoint=False))

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14, 7))
        batches = adata.obs[batch_key].astype('category')
        batch_categories = batches.cat.categories

        # Plot each scatter
        for i, (x, y, title) in enumerate([
            ('total_counts', 'n_genes_by_counts', 'n_genes vs total_counts'),
            ('total_counts', 'pct_counts_mt', 'pct_mt vs total_counts'),
        ]):
            for j, batch in enumerate(batch_categories):
                mask = batches == batch
                ax[i].scatter(adata.obs.loc[mask, x], adata.obs.loc[mask, y],
                             c=[colors[j]], s=5, alpha=0.6, label=batch if i == 0 else "")
            ax[i].set_xlabel(x)
            ax[i].set_ylabel(y)
            ax[i].set_title(title)

        # Add shared legend at bottom
        ncol = min(n_batches, 8)
        fig.legend(loc='lower center', ncol=ncol,
                   bbox_to_anchor=(0.5, -0.05), frameon=False, fontsize=8)
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.12 + 0.02 * (n_batches // 8))
        plt.savefig("{}/{}before_QC_scatter.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}before_QC_scatter.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()
    else:
        # Plot without grouping
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
        for i, feature in enumerate(features):
            sc.pl.violin(adata, keys=feature, stripplot=False, ax=axes.flatten()[i], show=False)
        plt.tight_layout()
        plt.savefig("{}/{}before_QC_violin.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}before_QC_violin.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14, 5))
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=ax[0], show=False)
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax[1], show=False)
        plt.tight_layout()
        plt.savefig("{}/{}before_QC_scatter.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}before_QC_scatter.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

    # Per-slice QC: derive n_genes / total_counts cutoffs from quantiles.
    # Rationale: stereo-seq slices differ in sequencing depth / tissue density,
    # so a global threshold either over-filters shallow slices or under-filters
    # deep ones. Using per-slice quantiles (default 1%/99%) keeps the bulk of
    # each slice while removing only the obvious outliers.
    print("Per-slice QC: n_genes / total_counts quantiles [{:.0%}, {:.0%}], mt<={}\n".format(
        quantile_low, quantile_high, max_mt))

    n_cells_before = adata.n_obs
    if batch_key in adata.obs.columns:
        # keep_mask starts all-False; each slice's keep decision is WRITTEN into
        # its own positions (NOT AND-accumulated, which would zero out other
        # slices' cells and leave 0 cells total).
        keep_mask = np.zeros(adata.n_obs, dtype=bool)
        threshold_table = []
        for sl in adata.obs[batch_key].astype('category').cat.categories:
            sl_mask = (adata.obs[batch_key] == sl).values
            n_sl = int(sl_mask.sum())
            if n_sl < 10:
                # Too few cells to estimate quantiles — keep all cells in this slice.
                keep_mask[sl_mask] = True
                threshold_table.append((sl, n_sl, n_sl, '-', '-', '-', '-'))
                continue
            ng_q = np.quantile(adata.obs.loc[sl_mask, 'n_genes_by_counts'], [quantile_low, quantile_high])
            tc_q = np.quantile(adata.obs.loc[sl_mask, 'total_counts'],    [quantile_low, quantile_high])
            ng_lo, ng_hi = int(ng_q[0]), int(ng_q[1])
            tc_lo, tc_hi = int(tc_q[0]), int(tc_q[1])
            sl_keep = (
                (adata.obs['n_genes_by_counts'] >= ng_lo) &
                (adata.obs['n_genes_by_counts'] <= ng_hi) &
                (adata.obs['total_counts']      >= tc_lo) &
                (adata.obs['total_counts']      <= tc_hi)
            ).values
            keep_mask[sl_mask] = sl_keep[sl_mask]
            n_kept = int(sl_keep[sl_mask].sum())
            threshold_table.append((sl, n_sl, n_kept, ng_lo, ng_hi, tc_lo, tc_hi))

        # Apply the global mt cutoff (mt% is meaningful per-cell, not per-slice).
        keep_mask &= (adata.obs['pct_counts_mt'] <= max_mt)

        # Print the per-slice threshold table
        print("Per-slice thresholds:")
        print("  {:<20s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}".format(
            "slice", "n_in", "n_keep", "ng_lo", "ng_hi", "tc_lo", "tc_hi"))
        for row in threshold_table:
            print("  {:<20s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}".format(*[str(x) for x in row]))
        print()
    else:
        # No batch_key — fall back to global quantiles (single-slice or unlabeled input).
        ng_q = np.quantile(adata.obs['n_genes_by_counts'], [quantile_low, quantile_high])
        tc_q = np.quantile(adata.obs['total_counts'],     [quantile_low, quantile_high])
        print("Global thresholds (no -BK): n_genes [{:.0f}, {:.0f}], total_counts [{:.0f}, {:.0f}], mt<={}\n".format(
            ng_q[0], ng_q[1], tc_q[0], tc_q[1], max_mt))
        keep_mask = (
            (adata.obs['n_genes_by_counts'] >= ng_q[0]) &
            (adata.obs['n_genes_by_counts'] <= ng_q[1]) &
            (adata.obs['total_counts']      >= tc_q[0]) &
            (adata.obs['total_counts']      <= tc_q[1]) &
            (adata.obs['pct_counts_mt']      <= max_mt)
        ).values

    print("Before filtering: {} cells x {} genes".format(adata.n_obs, adata.n_vars))
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[keep_mask].copy()
    print("After filtering : {} cells x {} genes\n".format(adata.n_obs, adata.n_vars))

    # Plot QC metrics after filtering
    if batch_key in adata.obs.columns:
        # Violin plots grouped by batch (1x3 grid, no stripplot)
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
        for ax, feature in zip(axes.flatten(), features):
            sc.pl.violin(adata, keys=feature, groupby=batch_key, stripplot=False,
                         rotation=45, ax=ax, show=False)
        plt.tight_layout()
        plt.savefig("{}/{}after_QC_violin.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}after_QC_violin.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

        # Scatter plots colored by batch (shared legend at bottom)
        n_batches = len(adata.obs[batch_key].astype('category').cat.categories)
        # Use enough colors for all batches
        if n_batches <= 20:
            colors = plt.cm.tab20(range(n_batches))
        elif n_batches <= 40:
            colors = np.vstack([plt.cm.tab20(range(20)), plt.cm.tab20b(range(n_batches - 20))])
        else:
            colors = plt.cm.hsv(np.linspace(0, 1, n_batches, endpoint=False))

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14, 7))
        batches = adata.obs[batch_key].astype('category')
        batch_categories = batches.cat.categories

        # Plot each scatter
        for i, (x, y, title) in enumerate([
            ('total_counts', 'n_genes_by_counts', 'n_genes vs total_counts'),
            ('total_counts', 'pct_counts_mt', 'pct_mt vs total_counts'),
        ]):
            for j, batch in enumerate(batch_categories):
                mask = batches == batch
                ax[i].scatter(adata.obs.loc[mask, x], adata.obs.loc[mask, y],
                             c=[colors[j]], s=5, alpha=0.6, label=batch if i == 0 else "")
            ax[i].set_xlabel(x)
            ax[i].set_ylabel(y)
            ax[i].set_title(title)

        # Add shared legend at bottom
        ncol = min(n_batches, 8)
        fig.legend(loc='lower center', ncol=ncol,
                   bbox_to_anchor=(0.5, -0.05), frameon=False, fontsize=8)
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.12 + 0.02 * (n_batches // 8))
        plt.savefig("{}/{}after_QC_scatter.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}after_QC_scatter.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()
    else:
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
        for i, feature in enumerate(features):
            sc.pl.violin(adata, keys=feature, stripplot=False, ax=axes.flatten()[i], show=False)
        plt.tight_layout()
        plt.savefig("{}/{}after_QC_violin.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}after_QC_violin.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14, 5))
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=ax[0], show=False)
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax[1], show=False)
        plt.tight_layout()
        plt.savefig("{}/{}after_QC_scatter.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}after_QC_scatter.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

    # ------------------------------------------------------------------ #
    # Step 3: Normalization (per-slice)                                     #
    # ------------------------------------------------------------------ #
    # Per-slice normalize_total + log1p, then reassemble into the merged
    # adata. Each slice's cells are normalized to the same library size
    # (target_sum=1e4) within that slice — prevents shallow slices from
    # being dragged down by deeper slices when normalize runs on the merged
    # matrix. counts / log1p layers are still kept on the merged adata for
    # downstream HVG / PCA / Harmony, identical to before.
    print(">>> Step 3: Normalization (per-slice)\n")

    adata.layers['counts'] = adata.X.copy()

    if batch_key in adata.obs.columns:
        n_slices = adata.obs[batch_key].astype('category').cat.categories.size
        print(f"Normalizing {n_slices} slices independently (target_sum=1e4 + log1p)\n")
        # Build a normalized matrix in slice order, then assign back.
        X_norm_chunks = []
        for sl in adata.obs[batch_key].astype('category').cat.categories:
            sl_mask = (adata.obs[batch_key] == sl).values
            sub = adata[sl_mask].copy()
            sc.pp.normalize_total(sub, target_sum=1e4)
            sc.pp.log1p(sub)
            X_norm_chunks.append(sub.X)
        # Concatenate in original obs order — chunk order matches the order
        # `cat.categories` walked, which may differ from obs order. Reorder
        # by reassigning via a single per-slice write into the merged X.
        X_new = adata.X.copy()
        for sl, chunk in zip(
            adata.obs[batch_key].astype('category').cat.categories,
            X_norm_chunks,
        ):
            sl_mask = (adata.obs[batch_key] == sl).values
            X_new[sl_mask] = chunk
        adata.X = X_new
    else:
        # No batch_key — single global normalize (same as before).
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    adata.layers['log1p'] = adata.X.copy()
    print("Per-slice normalized and log1p transformed\n")

    # ------------------------------------------------------------------ #
    # Step 4: Feature selection                                            #
    # ------------------------------------------------------------------ #
    print(">>> Step 4: Feature selection (HVGs)\n")

    adata.raw = adata
    if batch_key is None:
        sc.pp.highly_variable_genes(adata, layer="counts", n_top_genes=n_top_genes,
                                    flavor='seurat_v3', inplace=True)
    elif batch_key in adata.obs.columns:
        sc.pp.highly_variable_genes(adata, layer="counts", batch_key=batch_key,
                                    n_top_genes=n_top_genes, flavor='seurat_v3', inplace=True)
    else:
        print("Error: batch_key '{}' not found in adata.obs.columns".format(batch_key), file=sys.stderr)
        sys.exit(1)

    print("Found {} HVGs out of {} genes\n".format(adata.var['highly_variable'].sum(), adata.n_vars))

    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig("{}/{}HVGs.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}HVGs.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    # ------------------------------------------------------------------ #
    # Step 5: Scaling                                                      #
    # ------------------------------------------------------------------ #
    print(">>> Step 5: Scaling\n")
    # adata.layers["scaled"] = adata.X.toarray()
    # sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"], layer="scaled")
    # sc.pp.scale(adata, max_value=10, layer="scaled")
    sc.pp.scale(adata, max_value=10)
    print("Data scaled (max_value=10)\n")

    # ------------------------------------------------------------------ #
    # Step 6: PCA                                                          #
    # ------------------------------------------------------------------ #
    print(">>> Step 6: PCA\n")

    sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
    print("Computed 50 PCs\n")

    sc.pl.pca_loadings(adata, show=False)
    plt.savefig(f"{outpath}/{prefix}pca_loadings.png",
            dpi=300,
            bbox_inches='tight')
    plt.savefig(f"{outpath}/{prefix}pca_loadings.pdf",
            dpi=300,
            bbox_inches='tight')
    plt.close()

    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=False, show=False)
    plt.savefig("{}/{}pca_variance_ratio.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}pca_variance_ratio.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    if batch_key is not None:
        sc.pl.pca(adata, color=[batch_key, "pct_counts_mt"], dimensions=[(0, 1), (0, 1)],
                  ncols=2, size=2, wspace=0.4, show=False)
    else:
        sc.pl.pca(adata, color="pct_counts_mt", dimensions=[(0, 1)], size=2, show=False)
    plt.savefig("{}/{}pca.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}pca.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    # ------------------------------------------------------------------ #
    # Step 7: Batch correction (Harmony)                                   #
    # ------------------------------------------------------------------ #
    print(">>> Step 7: Batch correction\n")

    if batch_key is not None:
        try:
            import harmonypy
        except ImportError:
            print("Warning: harmonypy not installed, skipping batch correction\n")
        else:
            ho = harmonypy.run_harmony(adata.obsm['X_pca'], adata.obs, batch_key, max_iter_harmony=20)
            X_harmony = np.asarray(ho.Z_corr.T, dtype=np.float32)
            if X_harmony.shape[0] != adata.n_obs:
                X_harmony = np.asarray(ho.Z_corr, dtype=np.float32)
            adata.obsm['X_pca_harmony'] = X_harmony
            print("Harmony batch correction completed (batch_key={})\n".format(batch_key))
    else:
        print("No batch key specified, skipping batch correction\n")

    # ------------------------------------------------------------------ #
    # Step 8: Neighbor graph                                               #
    # ------------------------------------------------------------------ #
    print(">>> Step 8: Building neighbor graph\n")

    if 'X_pca_harmony' in adata.obsm:
        sc.pp.neighbors(adata, n_pcs=n_pcs, use_rep='X_pca_harmony')
        print("Using X_pca_harmony, n_pcs={}\n".format(n_pcs))
    else:
        sc.pp.neighbors(adata, n_pcs=n_pcs, use_rep='X_pca')
        print("Using X_pca, n_pcs={}\n".format(n_pcs))

    # ------------------------------------------------------------------ #
    # Step 9: Clustering                                                  #
    # ------------------------------------------------------------------ #
    print(">>> Step 9: Clustering (Leiden)\n")

    for res in [0.2, 0.5, 0.8, 1.0, 1.2, 1.5]:
        sc.tl.leiden(adata, flavor='igraph', resolution=res, key_added="leiden_{}".format(res))

    sc.tl.umap(adata)
    sc.pl.umap(adata, color=["leiden_0.2", "leiden_0.5", "leiden_0.8",
                             "leiden_1.0", "leiden_1.2", "leiden_1.5"],
               ncols=3, wspace=0.4, hspace=0.3, show=False)
    plt.savefig("{}/{}umap_clusters.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}umap_clusters.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    sc.tl.leiden(adata, flavor='igraph', resolution=resolution)
    n_clusters = len(adata.obs["leiden"].cat.categories)
    print("Final resolution={}, {} clusters detected\n".format(resolution, n_clusters))

    # ------------------------------------------------------------------ #
    # Step 10: UMAP visualization                                          #
    # ------------------------------------------------------------------ #
    print(">>> Step 10: UMAP visualization\n")

    sc.tl.umap(adata)
    if batch_key is not None:
        sc.pl.umap(adata, color=["leiden", batch_key], title=[f"Leiden (resolution: {resolution})", batch_key], wspace=0.4, show=False)
    else:
        sc.pl.umap(adata, color=["leiden"], title=f"Leiden (resolution: {resolution})", show=False) 
    plt.savefig("{}/{}umap_leiden.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}umap_leiden.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    # ------------------------------------------------------------------ #
    # Step 11: Marker genes                                                #
    # ------------------------------------------------------------------ #
    print(">>> Step 11: Marker gene identification\n")
    print("Method: Wilcoxon rank-sum test (CPU)\n")

    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, ncols=4, show=False)
    plt.savefig("{}/{}marker_genes.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}marker_genes.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    clusters = adata.obs["leiden"].cat.categories
    all_markers = pd.concat(
        [sc.get.rank_genes_groups_df(adata, group=cluster).assign(cluster=cluster)
         for cluster in clusters],
        ignore_index=True
    )
    all_markers.to_csv("{}/{}all_markers.csv".format(outpath, prefix), index=False)
    print("Marker genes saved to {}all_markers.csv\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 12: Save results                                                #
    # ------------------------------------------------------------------ #
    print(">>> Step 12: Saving results\n")

    output_file = "{}/{}preprocessed.h5ad".format(outpath, prefix)
    adata.write_h5ad(output_file, compression="gzip")
    print("Data saved to: {}".format(output_file))
    print("Final shape  : {} cells x {} genes".format(adata.n_obs, adata.n_vars))

    elapsed = time.time() - pipeline_start
    print("Total time   : {:.1f} s\n".format(elapsed))

    print("======================================================================")
    print("Pipeline completed successfully!")
    print("======================================================================\n")


if __name__ == "__main__":
    main()
