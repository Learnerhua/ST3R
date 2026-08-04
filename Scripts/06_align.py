#!/path/to/envs/spateo_env/bin/python
import sys, argparse, warnings, os, time
from datetime import datetime

# Suppress pkg_resources deprecation warning from spateo
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=UserWarning, module="spateo")

import torch
import spateo as st
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad

warnings.filterwarnings("ignore")


def plot_slices_2d(adata, slices, cluster_key, spatial_key, slices_key):
    """Plot slices in 2D grid."""
    labels = adata.obs[cluster_key].cat.categories.tolist()
    colors = adata.uns[f"{cluster_key}_colors"]
    palette = dict(zip(labels, colors))

    fig, _ = st.pl.slices_2d(
        slices=slices,
        slices_key=slices_key,
        label_key=cluster_key,
        spatial_key=spatial_key,
        height=5,
        ncols=6,
        center_coordinate=True,
        show_legend=True,
        legend_kwargs={
            "loc": "upper center",
            "bbox_to_anchor": (0.5, 0.1),
            "ncol": 5,
            "borderaxespad": 0,
            "frameon": False,
        },
        palette=palette,
        save_show_or_return="return",
    )
    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Spatial transcriptomics slice alignment using Spateo morpho alignment."
    )

    # Input/Output arguments (from notebook)
    parser.add_argument("-I", "--input", required=True, metavar="H5AD_FILE",
                        help="Path to input h5ad file")
    parser.add_argument("-O", "--outpath", type=str, default=".", metavar="OUTPUT_DIR",
                        help="Output directory (default: current directory)")
    parser.add_argument("-P", "--prefix", type=str, default="", metavar="PREFIX",
                        help="Prefix for output files (default: none)")

    # Key arguments (from notebook)
    parser.add_argument("-SPK", "--spatial_key", type=str, default="spatial", metavar="STR",
                        help="Key in adata.obsm for spatial coordinates (default: spatial)")
    parser.add_argument("-SLK", "--slices_key", type=str, default="slice_id", metavar="STR",
                        help="Key in adata.obs for slice identifier (default: slice_id)")
    parser.add_argument("-CK", "--cluster_key", type=str, default="squidpy_domains", metavar="STR",
                        help="Key in adata.obs for cluster annotation (default: squidpy_domains)")
    parser.add_argument("-KA", "--aligned_key", type=str, default="aligned_spatial", metavar="STR",
                        help="Key to add aligned spatial coordinates (default: aligned_spatial)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    input_file   = args.input
    outpath      = args.outpath
    prefix       = args.prefix
    spatial_key  = args.spatial_key
    slices_key   = args.slices_key
    cluster_key  = args.cluster_key
    aligned_key  = args.aligned_key

    # ------------------------------------------------------------------ #
    # Setup device                                                        #
    # ------------------------------------------------------------------ #
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # ------------------------------------------------------------------ #
    # Setup environment                                                    #
    # ------------------------------------------------------------------ #
    os.makedirs(outpath, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Print run information                                                #
    # ------------------------------------------------------------------ #
    print("\n====================== Spateo Slice Alignment =======================\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Working directory:", os.getcwd())
    print("Input file      :", input_file)
    print("Output directory:", outpath)
    print("Output prefix   :", prefix if prefix else "(none)")
    print("Device          :", device)
    print("Spateo version  :", st.__version__)
    print("\nKey parameters:")
    print("  spatial_key   :", spatial_key)
    print("  slices_key    :", slices_key)
    print("  cluster_key   :", cluster_key)
    print("  aligned_key   :", aligned_key)
    print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
    print("\n===================================================================\n")

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Step 1: Load data                                                    #
    # ------------------------------------------------------------------ #
    print(">>> Step 1: Loading data\n")
    adata = sc.read_h5ad(input_file)
    print("Data shape: {} cells x {} genes\n".format(adata.n_obs, adata.n_vars))

    # ------------------------------------------------------------------ #
    # Step 2: Set log1p matrix                                             #
    # ------------------------------------------------------------------ #
    print(">>> Step 2: Setting log1p matrix\n")
    adata.X = adata.layers["log1p"].copy()
    print("Using log1p transformed matrix\n")

    # ------------------------------------------------------------------ #
    # Step 3: Preprocess slices                                            #
    # ------------------------------------------------------------------ #
    print(">>> Step 3: Preprocessing slices\n")

    # Extract Z coordinate from slice_id
    def extract_number(s):
        parts = str(s).split('_')
        return int(parts[-1])

    adata.obs['Z'] = adata.obs[slices_key].apply(extract_number)

    # Sort slices by Z coordinate
    ordered_slices = (
        adata.obs
        .sort_values('Z')
        .drop_duplicates(slices_key)[slices_key]
        .tolist()
    )
    adata.obs[slices_key] = pd.Categorical(
        adata.obs[slices_key],
        categories=ordered_slices,
        ordered=True
    )
    slices = [adata[adata.obs[slices_key] == s].copy() for s in ordered_slices]
    print("Found {} slices\n".format(len(slices)))

    # ------------------------------------------------------------------ #
    # Step 4: Visualization before alignment                                #
    # ------------------------------------------------------------------ #
    print(">>> Step 4: Visualization before alignment\n")

    # First plot (legend position may be incorrect)
    fig1 = plot_slices_2d(
        adata,
        slices,
        cluster_key=cluster_key,
        spatial_key=spatial_key,
        slices_key=slices_key,
    )
    plt.close()

    # Second plot (legend position should be correct)
    fig1 = plot_slices_2d(
        adata,
        slices,
        cluster_key=cluster_key,
        spatial_key=spatial_key,
        slices_key=slices_key,
    )
    fig1.savefig(os.path.join(outpath, f"{prefix}squidpy_2Dslices.png"),
                dpi=300, bbox_inches='tight')
    fig1.savefig(os.path.join(outpath, f"{prefix}squidpy_2Dslices.pdf"),
                bbox_inches='tight')
    plt.close()
    print("Saved: {}squidpy_2Dslices.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 5: Alignment                                                    #
    # ------------------------------------------------------------------ #
    print(">>> Step 5: Running morpho alignment\n")

    transformation = st.align.morpho_align_transformation(
        models=slices,
        spatial_key=spatial_key,
        key_added=aligned_key,
        device=device,
        verbose=False,
        rep_layer='X_pca_harmony',
        rep_field='obsm',
        dissimilarity='cos',
    )
    print("Alignment completed\n")

    # ------------------------------------------------------------------ #
    # Step 6: Apply transformation                                         #
    # ------------------------------------------------------------------ #
    print(">>> Step 6: Applying transformation\n")

    aligned_slices = st.align.morpho_align_apply_transformation(
        models=slices,
        spatial_key=spatial_key,
        key_added=aligned_key,
        transformation=transformation,
    )
    print("Transformation applied\n")

    # ------------------------------------------------------------------ #
    # Step 7: Visualization after alignment                                 #
    # ------------------------------------------------------------------ #
    print(">>> Step 7: Visualization after alignment\n")

    fig2 = plot_slices_2d(
        adata,
        aligned_slices,
        cluster_key=cluster_key,
        spatial_key=aligned_key,
        slices_key=slices_key,
    )
    fig2.savefig(os.path.join(outpath, f"{prefix}aligned_2Dslices.png"),
                dpi=300, bbox_inches='tight')
    fig2.savefig(os.path.join(outpath, f"{prefix}aligned_2Dslices.pdf"),
                bbox_inches='tight')
    plt.close()
    print("Saved: {}aligned_2Dslices.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 8: Overlap visualization                                         #
    # ------------------------------------------------------------------ #
    print(">>> Step 8: Overlap visualization\n")

    slices_num = len(aligned_slices)
    mid = slices_num // 2

    unaligned_slices_overlap_plot = ad.concat(aligned_slices[mid-1:mid+2])
    unaligned_slices_overlap_plot.obsm['plot_spatial'] = unaligned_slices_overlap_plot.obsm['spatial']
    unaligned_slices_overlap_plot.obs['title'] = 'Unaligned'

    aligned_slices_overlap_plot = unaligned_slices_overlap_plot.copy()
    aligned_slices_overlap_plot.obsm['plot_spatial'] = aligned_slices_overlap_plot.obsm[aligned_key]
    aligned_slices_overlap_plot.obs['title'] = 'Spateo aligned'

    fig3 = plot_slices_2d(
        adata,
        [unaligned_slices_overlap_plot, aligned_slices_overlap_plot],
        cluster_key=cluster_key,
        spatial_key='plot_spatial',
        slices_key='title',
    )
    fig3.savefig(os.path.join(outpath, f"{prefix}aligned_2Dslices_overlap.png"),
                dpi=300, bbox_inches='tight')
    fig3.savefig(os.path.join(outpath, f"{prefix}aligned_2Dslices_overlap.pdf"),
                bbox_inches='tight')
    plt.close()
    print("Saved: {}aligned_2Dslices_overlap.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 9: Save results                                                 #
    # ------------------------------------------------------------------ #
    print(">>> Step 9: Saving results\n")

    aligned_adata = ad.concat(aligned_slices)
    output_file = os.path.join(outpath, f"{prefix}adata_aligned.h5ad")
    aligned_adata.write_h5ad(output_file, compression='gzip')
    print("Data saved to: {}".format(output_file))
    print("Final shape  : {} cells x {} genes".format(aligned_adata.n_obs, aligned_adata.n_vars))

    elapsed = time.time() - pipeline_start
    print("Total time   : {:.1f} s\n".format(elapsed))

    print("===================================================================")
    print("Alignment completed successfully!")
    print("===================================================================\n")


if __name__ == "__main__":
    main()