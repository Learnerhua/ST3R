#!/oldhome/ouyjh/miniforge3/envs/scanpy/bin/python
import sys, argparse, warnings, platform, os, math, time
from datetime import datetime
import squidpy as sq
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=FutureWarning, module="scanpy")


def main():
    parser = argparse.ArgumentParser(
        description="Spatial transcriptomics analysis using Squidpy for domain detection."
    )

    parser.add_argument(
        "-I", "--input",
        required=True,
        metavar="H5AD_FILE",
        help="Path to the input H5AD file (preprocessed)"
    )

    parser.add_argument(
        "-O", "--outpath",
        type=str,
        default="Output",
        metavar="OUTPUT_DIR",
        help="Output directory for results (default: Output)"
    )

    parser.add_argument(
        "-P", "--prefix",
        type=str,
        default="",
        metavar="PREFIX",
        help="Prefix for output files (default: none)"
    )

    parser.add_argument(
        "-LK", "--library_key",
        type=str,
        default="slice_id",
        metavar="LIBRARY_KEY",
        help="Key in adata.obs for library/slice identifier (default: slice_id)"
    )

    parser.add_argument(
        "-A", "--alpha",
        type=float,
        default=0.2,
        metavar="ALPHA",
        help="Weight for spatial graph in joint graph: (1-alpha)*gene_graph + alpha*spatial_graph (default: 0.2)"
    )

    parser.add_argument(
        "-R", "--resolution",
        type=float,
        default=0.7,
        metavar="RESOLUTION",
        help="Resolution parameter for Leiden clustering (default: 0.7)"
    )

    parser.add_argument(
        "-NC", "--n_cols",
        type=int,
        default=6,
        metavar="N_COLS",
        help="Number of columns in the output grid plot (default: 6)"
    )

    parser.add_argument(
        "-WS", "--wspace",
        type=float,
        default=0.4,
        metavar="WSPACE",
        help="Width space between subplots (default: 0.4)"
    )

    parser.add_argument(
        "-HS", "--hspace",
        type=float,
        default=0.2,
        metavar="HSPACE",
        help="Height space between subplots (default: 0.2)"
    )

    parser.add_argument(
        "-SS", "--spot_size",
        type=float,
        default=15,
        metavar="SPOT_SIZE",
        help="Spot size for spatial visualization (default: 15)"
    )

    parser.add_argument(
        "-DPI", "--dpi",
        type=int,
        default=300,
        metavar="DPI",
        help="DPI for saving PNG images (default: 300)"
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    input_file   = args.input
    outpath      = args.outpath
    prefix       = args.prefix
    library_key  = args.library_key
    sq_alpha     = args.alpha
    resolution   = args.resolution
    n_cols       = args.n_cols
    wspace       = args.wspace
    hspace       = args.hspace
    spot_size    = args.spot_size
    dpi          = args.dpi

    # ------------------------------------------------------------------ #
    # Print run information                                                #
    # ------------------------------------------------------------------ #
    print("\n========================= Squidpy Spatial Analysis =========================\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Operating system:", platform.system(), platform.release())
    print("Platform        :", platform.platform())
    print("Working directory:", os.getcwd())
    print("Input H5AD file :", input_file)
    print("Output directory:", outpath)
    print("Output prefix   :", prefix if prefix else "(none)")
    print("library_key     :", library_key)
    print("alpha           :", sq_alpha)
    print("resolution      :", resolution)
    print("n_cols          :", n_cols)
    print("wspace          :", wspace)
    print("hspace          :", hspace)
    print("spot_size       :", spot_size)
    print("dpi             :", dpi)
    print("\nCopyright (c) 2026 OYJH. All Rights Reserved.")
    print("\n=============================================================================\n")

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Validate input file                                                  #
    # ------------------------------------------------------------------ #
    if not os.path.isfile(input_file):
        print(f"Error: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # Create output directory if needed                                    #
    # ------------------------------------------------------------------ #
    os.makedirs(outpath, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Step 1: Read H5AD file                                               #
    # ------------------------------------------------------------------ #
    print(f">>> Reading H5AD file: {input_file}\n")

    adata = sc.read_h5ad(input_file)
    print(adata)

    # ------------------------------------------------------------------ #
    # Step 2: Build spatial neighbors graph                                #
    # ------------------------------------------------------------------ #
    print("\n>>> Building spatial neighbors graph...\n")

    nn_graph_genes = adata.obsp["connectivities"]
    sq.gr.spatial_neighbors(adata, library_key=library_key)
    nn_graph_space = adata.obsp["spatial_connectivities"]

    # ------------------------------------------------------------------ #
    # Step 3: Construct joint graph                                        #
    # ------------------------------------------------------------------ #
    print(f"\n>>> Constructing joint graph with alpha={sq_alpha}...\n")

    joint_graph = (1 - sq_alpha) * nn_graph_genes + sq_alpha * nn_graph_space

    # ------------------------------------------------------------------ #
    # Step 4: Leiden clustering                                            #
    # ------------------------------------------------------------------ #
    print(f"\n>>> Running Leiden clustering with resolution={resolution}...\n")

    sc.tl.leiden(
        adata,
        adjacency=joint_graph,
        key_added="squidpy_domains",
        flavor="igraph",
        n_iterations=2,
        resolution=resolution
    )

    # ------------------------------------------------------------------ #
    # Step 4.5: UMAP visualization                                         #
    # ------------------------------------------------------------------ #
    sc.pl.umap(adata, color=["leiden", "squidpy_domains"], legend_loc='on data',
               title=["leiden", "squidpy"], show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, f"{prefix}leiden_squidpy.png"), dpi=dpi, bbox_inches="tight")
    plt.savefig(os.path.join(outpath, f"{prefix}leiden_squidpy.pdf"), bbox_inches="tight")
    plt.close()

    # ------------------------------------------------------------------ #
    # Step 5: Visualization                                                #
    # ------------------------------------------------------------------ #
    print("\n>>> Generating spatial visualization...\n")

    slices = adata.obs[library_key].cat.categories
    n_rows = math.ceil(len(slices) / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows), gridspec_kw={'wspace': wspace, 'hspace': hspace})
    axes = axes.flatten()

    for i, s in enumerate(slices):
        adata_tmp = adata[adata.obs[library_key] == s].copy()
        sc.pl.spatial(
            adata_tmp,
            basis="spatial",
            color="squidpy_domains",
            spot_size=spot_size,
            title=s,
            ax=axes[i],
            show=False
        )

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()

    png_path = os.path.join(outpath, f"{prefix}squidpy_domains_grid.png")
    pdf_path = os.path.join(outpath, f"{prefix}squidpy_domains_grid.pdf")
    plt.savefig(png_path, dpi=dpi, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")

    # Leiden grid plot
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows), gridspec_kw={'wspace': wspace, 'hspace': hspace})
    axes = axes.flatten()

    for i, s in enumerate(slices):
        adata_tmp = adata[adata.obs[library_key] == s].copy()
        sc.pl.spatial(
            adata_tmp,
            basis="spatial",
            color="leiden",
            spot_size=spot_size,
            title=s,
            ax=axes[i],
            show=False
        )

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()

    png_path = os.path.join(outpath, f"{prefix}leiden_grid.png")
    pdf_path = os.path.join(outpath, f"{prefix}leiden_grid.pdf")
    plt.savefig(png_path, dpi=dpi, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")

    # ------------------------------------------------------------------ #
    # Step 6: Save results                                                 #
    # ------------------------------------------------------------------ #
    h5ad_path = os.path.join(outpath, f"{prefix}squidpy.h5ad")
    adata.write_h5ad(h5ad_path, compression="gzip")

    print(f"\n>>> Results saved to: {h5ad_path}\n")
    print("Final shape  : {} cells x {} genes".format(adata.n_obs, adata.n_vars))

    elapsed = time.time() - pipeline_start
    print("Total time   : {:.1f} s\n".format(elapsed))

    print("=============================================================================")
    print("Analysis completed successfully!")
    print("=============================================================================\n")


if __name__ == "__main__":
    main()
