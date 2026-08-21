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
        description="Single-cell RNA-seq analysis pipeline (downstream): "
                    "HVG selection -> Scaling -> PCA -> Batch correction -> "
                    "Neighbor graph -> Leiden clustering -> UMAP -> Marker gene identification. "
                    "Per-sample QC and normalization are performed upstream by 02_concat.py."
    )

    # Input/Output arguments
    parser.add_argument("-I", "--input", required=True, metavar="H5AD_FILE",
                        help="Path to input h5ad file (concatenated log1p-normalized adata from 02_concat.py)")
    parser.add_argument("-O", "--outpath", type=str, default=".", metavar="OUTPUT_DIR",
                        help="Output directory (default: current directory)")
    parser.add_argument("-P", "--prefix", type=str, default="", metavar="PREFIX",
                        help="Prefix for output files (default: none)")

    # Batch correction arguments
    parser.add_argument("-BK", "--batch_key", type=str, default="batch", metavar="STR",
                        help="Batch key in adata.obs for batch correction (default: batch)")

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

    input_file   = args.input
    outpath      = args.outpath
    prefix       = args.prefix
    batch_key    = args.batch_key
    n_pcs        = args.n_pcs
    n_top_genes  = args.n_top_genes
    resolution   = args.resolution
    dpi          = args.dpi

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
    print("\n==================== scRNA-seq Analysis Pipeline (Downstream) ====================\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Working directory:", os.getcwd())
    print("Input file      :", input_file)
    print("Output directory:", outpath)
    print("Output prefix   :", prefix)
    print("Batch key       :", batch_key)
    print("Clustering      : n_pcs={}, n_top_genes={}, resolution={}".format(
          n_pcs, n_top_genes, resolution))
    print("\nNote: Per-sample QC and normalization are performed upstream by 02_concat.py.")
    print("      This script expects adata.X to already be log1p-normalized.")
    print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
    print("\n======================================================================\n")

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Step 1: Load data                                                    #
    # ------------------------------------------------------------------ #
    print(">>> Step 1: Loading data\n")
    adata = sc.read_h5ad(input_file)
    print("Data shape: {} cells x {} genes".format(adata.n_obs, adata.n_vars))
    print("X dtype   : {} (expect: float = log1p-normalized)".format(adata.X.dtype))

    # ------------------------------------------------------------------ #
    # Step 2: Feature selection (HVGs)                                     #
    # ------------------------------------------------------------------ #
    print("\n>>> Step 2: Feature selection (HVGs)\n")

    adata.raw = adata
    if batch_key in adata.obs.columns:
        sc.pp.highly_variable_genes(adata, layer="counts" if "counts" in adata.layers else None,
                                    batch_key=batch_key,
                                    n_top_genes=n_top_genes, flavor='seurat_v3', inplace=True)
    else:
        sc.pp.highly_variable_genes(adata,
                                    n_top_genes=n_top_genes, flavor='seurat_v3', inplace=True)

    print("Found {} HVGs out of {} genes\n".format(adata.var['highly_variable'].sum(), adata.n_vars))

    # Plot HVG scatter (filename: {prefix}HVGs.{png,pdf}) — required by Report_config
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig("{}/{}HVGs.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}HVGs.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    # ------------------------------------------------------------------ #
    # Step 3: Scaling                                                      #
    # ------------------------------------------------------------------ #
    print(">>> Step 3: Scaling\n")
    # adata.layers["scaled"] = adata.X.toarray()
    # sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"], layer="scaled")
    # sc.pp.scale(adata, max_value=10, layer="scaled")
    sc.pp.scale(adata, max_value=10)
    print("Data scaled (max_value=10)\n")

    # ------------------------------------------------------------------ #
    # Step 4: PCA                                                          #
    # ------------------------------------------------------------------ #
    print(">>> Step 4: PCA\n")

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

    if batch_key in adata.obs.columns:
        sc.pl.pca(adata, color=[batch_key, "pct_counts_mt"], dimensions=[(0, 1), (0, 1)],
                  ncols=2, size=2, wspace=0.4, show=False)
    else:
        sc.pl.pca(adata, color="pct_counts_mt", dimensions=[(0, 1)], size=2, show=False)
    plt.savefig("{}/{}pca.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}pca.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    # ------------------------------------------------------------------ #
    # Step 5: Batch correction (Harmony)                                   #
    # ------------------------------------------------------------------ #
    print(">>> Step 5: Batch correction\n")

    if batch_key in adata.obs.columns:
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
    # Step 6: Neighbor graph                                               #
    # ------------------------------------------------------------------ #
    print(">>> Step 6: Building neighbor graph\n")

    if 'X_pca_harmony' in adata.obsm:
        sc.pp.neighbors(adata, n_pcs=n_pcs, use_rep='X_pca_harmony')
        print("Using X_pca_harmony, n_pcs={}\n".format(n_pcs))
    else:
        sc.pp.neighbors(adata, n_pcs=n_pcs, use_rep='X_pca')
        print("Using X_pca, n_pcs={}\n".format(n_pcs))

    # ------------------------------------------------------------------ #
    # Step 7: Clustering (Leiden)                                          #
    # ------------------------------------------------------------------ #
    print(">>> Step 7: Clustering (Leiden)\n")

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
    # Step 8: UMAP visualization                                           #
    # ------------------------------------------------------------------ #
    print(">>> Step 8: UMAP visualization\n")

    sc.tl.umap(adata)
    if batch_key in adata.obs.columns:
        sc.pl.umap(adata, color=["leiden", batch_key], title=[f"Leiden (resolution: {resolution})", batch_key], wspace=0.4, show=False)
    else:
        sc.pl.umap(adata, color=["leiden"], title=f"Leiden (resolution: {resolution})", show=False)
    plt.savefig("{}/{}umap_leiden.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}umap_leiden.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    # ------------------------------------------------------------------ #
    # Step 9: Marker genes                                                 #
    # ------------------------------------------------------------------ #
    print(">>> Step 9: Marker gene identification\n")
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
    # Step 10: Save results                                                #
    # ------------------------------------------------------------------ #
    print(">>> Step 10: Saving results\n")

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