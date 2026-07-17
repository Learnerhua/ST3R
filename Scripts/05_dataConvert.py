#!/oldhome/ouyjh/miniforge3/envs/scanpy/bin/python
import sys, argparse, os, time
from datetime import datetime
import scanpy as sc


def main():
    parser = argparse.ArgumentParser(
        description="Clean h5ad file by removing specific uns and obsp keys"
    )

    # Input/Output arguments
    parser.add_argument("-I", "--input", required=True, metavar="H5AD_FILE",
                        help="Path to input h5ad file")
    parser.add_argument("-O", "--outpath", type=str, default=".", metavar="OUTPUT_DIR",
                        help="Output directory (default: current directory)")
    parser.add_argument("-P", "--prefix", type=str, default="", metavar="PREFIX",
                        help="Prefix for output file (default: none)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    input_file = args.input
    outpath = args.outpath
    prefix = args.prefix

    # Set output file path
    output_file = os.path.join(outpath, f"{prefix}compatible.h5ad")

    # Create output directory if it doesn't exist
    os.makedirs(outpath, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Print run information                                                #
    # ------------------------------------------------------------------ #
    print("\n==================== h5ad Data Cleaning ====================\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Working directory:", os.getcwd())
    print("Input file      :", input_file)
    print("Output directory:", outpath)
    print("Output prefix   :", prefix if prefix else "(none)")
    print("\nCopyright (c) 2026 OYJH. All Rights Reserved.")
    print("\n======================================================================\n")

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Step 1: Load data                                                    #
    # ------------------------------------------------------------------ #
    print(">>> Step 1: Loading data\n")
    adata = sc.read_h5ad(input_file)
    print("Data shape: {} cells x {} genes\n".format(adata.n_obs, adata.n_vars))

    # ------------------------------------------------------------------ #
    # Step 2: Remove specific keys                                         #
    # ------------------------------------------------------------------ #
    print(">>> Step 2: Removing specific keys\n")

    # Keys to remove from .uns
    uns_keys_to_remove = ['log1p', 'rank_genes_groups', 'spatial_neighbors', 'squidpy_domains']

    # Keys to remove from .obsp
    obsp_keys_to_remove = ['spatial_connectivities', 'spatial_distances']

    # Remove uns keys
    for key in uns_keys_to_remove:
        if key in adata.uns.keys():
            print("Removing .uns['{}']".format(key))
            del adata.uns[key]
        else:
            print("Warning: .uns['{}'] not found, skipping".format(key))

    # Remove obsp keys
    for key in obsp_keys_to_remove:
        if key in adata.obsp.keys():
            print("Removing .obsp['{}']".format(key))
            del adata.obsp[key]
        else:
            print("Warning: .obsp['{}'] not found, skipping".format(key))

    print()

    # ------------------------------------------------------------------ #
    # Step 3: Save cleaned data                                            #
    # ------------------------------------------------------------------ #
    print(">>> Step 3: Saving cleaned data\n")

    adata.write_h5ad(output_file, compression='gzip')
    print("Data saved to: {}".format(output_file))
    print("Final shape : {} cells x {} genes".format(adata.n_obs, adata.n_vars))

    elapsed = time.time() - pipeline_start
    print("Total time  : {:.1f} s\n".format(elapsed))

    print("======================================================================")
    print("Data cleaning completed successfully!")
    print("======================================================================\n")


if __name__ == "__main__":
    main()
