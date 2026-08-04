#!/path/to/envs/stereopy/bin/python
import sys, argparse, warnings, platform, os, re, time
from datetime import datetime
import pandas as pd
import stereo as st


def process_single_file(input_file, outpath, prefix, bin_type, bin_size, threads, flavor, sample_id):
    """Process a single GEF file (original behavior)."""

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Validate input file                                                  #
    # ------------------------------------------------------------------ #
    if not input_file:
        print("Error: Input file is required in single file mode.", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(input_file):
        print(f"Error: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # Derive output filename from input stem                               #
    # bins mode:      abcd.gef + -BS 50  -> abcd_bin50.h5ad               #
    # cell_bins mode: abcd.gef           -> abcd.h5ad                     #
    # ------------------------------------------------------------------ #
    input_stem = os.path.splitext(os.path.basename(input_file))[0]
    if bin_type == "bins":
        output_name = f"{prefix}{input_stem}_bin{bin_size}.h5ad"
    else:
        output_name = f"{prefix}{input_stem}.h5ad"
    output_path = os.path.join(outpath, output_name)

    # ------------------------------------------------------------------ #
    # Print run information                                                #
    # ------------------------------------------------------------------ #
    print("\n========================= GEF to H5AD Conversion =========================\n")
    print("Current time   :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Operating system:", platform.system(), platform.release())
    print("Platform       :", platform.platform())
    print("Working directory:", os.getcwd())
    print("Input GEF file :", input_file)
    print("Output H5AD    :", output_path)
    print("bin_type       :", bin_type)
    if bin_type == "bins":
        print("bin_size       :", bin_size)
    print("threads        :", threads)
    print("flavor         :", flavor)
    print("sample_id      :", sample_id)
    print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
    print("\n===========================================================================\n")

    # ------------------------------------------------------------------ #
    # Create output directory if needed                                    #
    # ------------------------------------------------------------------ #
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # ------------------------------------------------------------------ #
    # Step 1: Read GEF file                                                #
    # ------------------------------------------------------------------ #
    print(f">>> Reading GEF file: {input_file}\n")

    read_kwargs = dict(
        file_path=input_file,
        bin_type=bin_type,
        is_sparse=True,
        gene_name_index=True,
        num_threads=threads,
    )
    if bin_type == "bins":
        read_kwargs["bin_size"] = bin_size

    data = st.io.read_gef(**read_kwargs)

    print("\nStereoExpData object:\n")
    print(data)

    # ------------------------------------------------------------------ #
    # Step 2: Convert to AnnData and save as H5AD                         #
    # ------------------------------------------------------------------ #
    print(f"\n>>> Converting to AnnData and saving to: {output_path}\n")

    st.io.stereo_to_anndata(
        data=data,
        flavor=flavor,
        sample_id=sample_id,
        reindex=False,
        output=output_path,
        compression="gzip",
    )

    elapsed = time.time() - pipeline_start
    print(f"\nFile: {output_name} saved successfully !")
    print("Total time   : {:.1f} s\n".format(elapsed))


def process_config_mode(config_df, outpath, prefix, bin_type, bin_size, threads, flavor):
    """Process multiple GEF files from config file."""

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Create output directory if needed                                    #
    # ------------------------------------------------------------------ #
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # ------------------------------------------------------------------ #
    # Process each GEF file                                               #
    # ------------------------------------------------------------------ #
    success_count = 0
    error_count = 0

    for idx, row in config_df.iterrows():
        gef_path = row["gef_path"]
        sample_id = row["sample_id"]
        slice_id = row["slice_id"]

        # Derive output filename
        input_stem = os.path.splitext(os.path.basename(gef_path))[0]
        if bin_type == "bins":
            output_name = f"{prefix}{input_stem}_bin{bin_size}.h5ad"
        else:
            output_name = f"{prefix}{input_stem}.h5ad"
        output_path = os.path.join(outpath, output_name)

        print(f"\n[{idx+1}/{len(config_df)}] Processing: {gef_path}")
        print(f"    sample_id: {sample_id}, slice_id: {slice_id}")
        if bin_type == "bins":
            print(f"    bin_type: {bin_type}, bin_size: {bin_size}")
        else:
            print(f"    bin_type: {bin_type}")
        print(f"    output: {output_name}")

        try:
            # ------------------------------------------------------------------ #
            # Step 1: Read GEF file                                                #
            # ------------------------------------------------------------------ #
            read_kwargs = dict(
                file_path=gef_path,
                bin_type=bin_type,
                is_sparse=True,
                gene_name_index=True,
                num_threads=threads,
            )
            if bin_type == "bins" and bin_size:
                read_kwargs["bin_size"] = bin_size

            data = st.io.read_gef(**read_kwargs)

            # ------------------------------------------------------------------ #
            # Step 2: Convert to AnnData                                          #
            # ------------------------------------------------------------------ #
            adata = st.io.stereo_to_anndata(
                data=data,
                flavor=flavor,
                sample_id=sample_id,
                reindex=False,
                output=None,
            )

            # Add slice_id to obs
            adata.obs["slice_id"] = slice_id

            # ------------------------------------------------------------------ #
            # Step 3: Save as H5AD                                                #
            # ------------------------------------------------------------------ #
            adata.write_h5ad(output_path, compression="gzip")

            print(f"    [SUCCESS] Saved: {output_name}")
            success_count += 1

        except Exception as e:
            print(f"    [ERROR] {str(e)}")
            error_count += 1

    # ------------------------------------------------------------------ #
    # Summary                                                              #
    # ------------------------------------------------------------------ #
    elapsed = time.time() - pipeline_start
    print("\n===========================================================================")
    print(f"Batch processing completed!")
    print(f"Success: {success_count}/{len(config_df)}")
    print(f"Errors: {error_count}/{len(config_df)}")
    print("Total time   : {:.1f} s".format(elapsed))
    print("===========================================================================\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert GEF file to H5AD format using Stereo-seq pipeline."
    )

    parser.add_argument(
        "-I", "--input",
        type=str,
        default=None,
        metavar="GEF_FILE",
        help="Path to the input GEF file (mutually exclusive with -C)"
    )

    parser.add_argument(
        "-C", "--config",
        type=str,
        default=None,
        metavar="CONFIG_FILE",
        help="Config file for batch processing. TSV format with header row containing columns: "
             "gef_path (path to GEF file), sample_id (sample identifier), slice_id (slice identifier). "
             "When -C is specified, -I/-S are ignored and -BT/-BS apply to all files. "
    )

    parser.add_argument(
        "-O", "--outpath",
        type=str,
        default=".",
        metavar="OUTPUT_DIR",
        help="Output directory for the H5AD file (default: current directory)"
    )

    parser.add_argument(
        "-BT", "--bin_type",
        type=str,
        default="bins",
        choices=["bins", "cell_bins"],
        metavar="BIN_TYPE",
        help="Bin type: 'bins' for fixed-grid bins, 'cell_bins' for cell-segmented data (default: bins)"
    )

    parser.add_argument(
        "-BS", "--bin_size",
        type=int,
        default=50,
        metavar="BIN_SIZE",
        help="Bin size, only effective when --bin_type is 'bins' (default: 50)"
    )

    parser.add_argument(
        "-P", "--prefix",
        type=str,
        default="",
        metavar="PREFIX",
        help="Prefix for the output H5AD filename (default: none)"
    )

    parser.add_argument(
        "-T", "--threads",
        type=int,
        default=-1,
        metavar="THREADS",
        help="Number of threads to use; -1 means all available CPU cores (default: -1)"
    )

    parser.add_argument(
        "-F", "--flavor",
        type=str,
        default="scanpy",
        choices=["scanpy", "seurat"],
        metavar="FLAVOR",
        help="AnnData flavor: 'scanpy' for Python ecosystem, 'seurat' for R Seurat compatibility (default: scanpy)"
    )

    parser.add_argument(
        "-S", "--sample_id",
        type=str,
        default="sample",
        metavar="SAMPLE_ID",
        help="Sample ID written into adata.obs['orig.ident'] (default: sample)"
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    input_file  = args.input
    outpath     = args.outpath
    prefix      = args.prefix
    bin_type    = args.bin_type
    bin_size    = args.bin_size
    threads     = args.threads
    flavor      = args.flavor
    sample_id   = args.sample_id
    config_file = args.config

    # ------------------------------------------------------------------ #
    # Parameter validation                                                #
    # ------------------------------------------------------------------ #
    if not config_file and not input_file:
        print("Error: Either -I or -C must be specified.", file=sys.stderr)
        sys.exit(1)

    if config_file and input_file:
        print("Error: Cannot use both -I and -C options at the same time.", file=sys.stderr)
        sys.exit(1)

    if config_file:
        # Config file mode
        if not os.path.isfile(config_file):
            print(f"Error: Config file not found: {config_file}", file=sys.stderr)
            sys.exit(1)

        # Print run information first
        print("\n========================= Batch GEF to H5AD Conversion =========================\n")
        print("Current time   :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print("Operating system:", platform.system(), platform.release())
        print("Platform       :", platform.platform())
        print("Working directory:", os.getcwd())
        print("Config file    :", config_file)
        print("Output directory:", outpath)
        print("bin_type       :", bin_type)
        if bin_type == "bins":
            print("bin_size       :", bin_size)
        print("threads        :", threads)
        print("flavor         :", flavor)
        print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
        print("\n===========================================================================\n")

        # Then load config file
        print(f"\n>>> Loading config file: {config_file}")
        config_df = pd.read_csv(config_file, sep="\t", header=0)

        # Validate columns
        required_cols = ["gef_path", "sample_id", "slice_id"]
        if not all(col in config_df.columns for col in required_cols):
            print(f"Error: Config file must contain columns: {', '.join(required_cols)}", file=sys.stderr)
            print(f"Found columns: {', '.join(config_df.columns)}", file=sys.stderr)
            sys.exit(1)

        print(f"Found {len(config_df)} entries in config file")
        print("Total files    :", len(config_df))
        print(config_df.head())

        # Validate all GEF files exist
        missing_files = []
        for idx, row in config_df.iterrows():
            if not os.path.isfile(row["gef_path"]):
                missing_files.append(row["gef_path"])

        if missing_files:
            print(f"Error: {len(missing_files)} GEF file(s) not found:", file=sys.stderr)
            for f in missing_files:
                print(f"  - {f}", file=sys.stderr)
            sys.exit(1)

        process_config_mode(config_df, outpath, prefix, bin_type, bin_size, threads, flavor)
    else:
        # Single file mode (original behavior)
        process_single_file(
            input_file, outpath, prefix, bin_type, bin_size,
            threads, flavor, sample_id
        )


if __name__ == "__main__":
    main()
