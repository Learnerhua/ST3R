#!/path/to/envs/scanpy/bin/python
import sys, argparse, warnings, platform, os, glob, time
from datetime import datetime
import anndata as ad


def concat_h5ad(input_dir, outpath, prefix, join_method, merge_method, uns_merge_method, label_name):
    """Concatenate multiple H5AD files from a directory."""

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
    # Load all H5AD files                                                  #
    # ------------------------------------------------------------------ #
    print("\n>>> Loading H5AD files...")
    adata_list = []
    sample_names = []

    for i, h5ad_file in enumerate(h5ad_files):
        try:
            adata = ad.read_h5ad(h5ad_file)
            adata_list.append(adata)

            # Get sample name from slice_id if available, otherwise from filename
            if "slice_id" in adata.obs.columns:
                sample_name = adata.obs["slice_id"].iloc[0]
            elif "orig.ident" in adata.obs.columns:
                sample_name = adata.obs["orig.ident"].iloc[0]
            else:
                sample_name = os.path.splitext(os.path.basename(h5ad_file))[0]

            sample_names.append(sample_name)
            print(f"  [{i+1}/{len(h5ad_files)}] {os.path.basename(h5ad_file)}: "
                  f"{adata.n_obs} cells, {adata.n_vars} genes, sample: {sample_name}")

        except Exception as e:
            print(f"  [ERROR] Failed to load {h5ad_file}: {str(e)}", file=sys.stderr)
            sys.exit(1)

    # ------------------------------------------------------------------ #
    # Concatenate H5AD files                                               #
    # ------------------------------------------------------------------ #
    print(f"\n>>> Concatenating {len(adata_list)} files...")
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
        description="Concatenate multiple H5AD files into a single file."
    )

    parser.add_argument(
        "-I", "--input",
        type=str,
        required=True,
        metavar="INPUT_DIR",
        help="Directory containing H5AD files to concatenate"
    )

    parser.add_argument(
        "-O", "--outpath",
        type=str,
        default=".",
        metavar="OUTPUT_DIR",
        help="Output directory for the concatenated H5AD file (default: current directory)"
    )

    parser.add_argument(
        "-P", "--prefix",
        type=str,
        default="",
        metavar="PREFIX",
        help="Prefix for the output H5AD filename (default: none)"
    )

    parser.add_argument(
        "--join",
        type=str,
        default="outer",
        choices=["outer", "inner"],
        metavar="JOIN",
        help="Join method: 'outer' keeps all genes, 'inner' keeps only common genes (default: outer)"
    )

    parser.add_argument(
        "--merge",
        type=str,
        default="same",
        choices=["same", "unique", "first", "outer", "none"],
        metavar="MERGE",
        help="Merge method for obs/var: 'same' keeps common columns, 'unique' keeps all unique columns (default: same)"
    )

    parser.add_argument(
        "--uns-merge",
        type=str,
        default="unique",
        choices=["same", "unique", "first", "outer", "none"],
        metavar="UNS_MERGE",
        help="Merge method for uns dict keys: 'unique' keeps all unique keys (default: unique)"
    )

    parser.add_argument(
        "--label",
        type=str,
        default="batch",
        metavar="LABEL",
        help="Column name for batch labels in obs (default: batch)"
    )

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

    # ------------------------------------------------------------------ #
    # Print run information                                                #
    # ------------------------------------------------------------------ #
    print("\n=========================== H5AD Concatenation ============================\n")
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
    print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
    print("\n===========================================================================\n")

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
    # Concatenate H5AD files                                               #
    # ------------------------------------------------------------------ #
    concat_h5ad(input_dir, outpath, prefix, join_method, merge_method, uns_merge_method, label_name)


if __name__ == "__main__":
    main()
