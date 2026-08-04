#!/path/to/envs/spateo_env/bin/python
import sys, argparse, warnings, os, time
from datetime import datetime

# Suppress pkg_resources deprecation warning from spateo
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=UserWarning, module="spateo")

import numpy as np
import spateo as st
import scanpy as sc
import pandas as pd
import scipy.sparse as sp
import pyvista as pv

# Start virtual display for off-screen rendering on servers
pv.start_xvfb()

warnings.filterwarnings("ignore")


def main():
    parser = argparse.ArgumentParser(
        description="Gaussian Process interpolation for gene expression in 3D spatial transcriptomics data."
    )

    # Input/Output arguments
    parser.add_argument("-AD", "--aligned_data", required=True, metavar="H5AD_FILE",
                        help="Path to aligned h5ad file from TDR step")
    parser.add_argument("-PC", "--pc_model", required=True, metavar="VTK_FILE",
                        help="Path to point cloud model VTK file")
    parser.add_argument("-MS", "--mesh_model", required=True, metavar="VTK_FILE",
                        help="Path to mesh model VTK file")
    parser.add_argument("-VX", "--voxel_model", required=True, metavar="VTK_FILE",
                        help="Path to voxel model VTK file")
    parser.add_argument("-GL", "--glm_results", required=True, metavar="CSV_FILE",
                        help="Path to GLM DE results CSV file")
    parser.add_argument("-O", "--outpath", type=str, default=".", metavar="OUTPUT_DIR",
                        help="Output directory (default: current directory)")
    parser.add_argument("-P", "--prefix", type=str, default="", metavar="PREFIX",
                        help="Prefix for output files (default: none)")

    # Key arguments
    parser.add_argument("-TK", "--tdr_key", type=str, default="aligned_spatial_3D", metavar="STR",
                        help="Key for 3D spatial coordinates (default: aligned_spatial_3D)")
    parser.add_argument("-NG", "--num_genes", type=int, default=3, metavar="INT",
                        help="Number of top genes from GLM results to interpolate (default: 3)")

    # Slice parameters
    parser.add_argument("-NS", "--n_slices", type=int, default=15, metavar="INT",
                        help="Number of slices for 3D slice visualization (default: 15)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    aligned_data  = args.aligned_data
    pc_model      = args.pc_model
    mesh_model    = args.mesh_model
    voxel_model   = args.voxel_model
    glm_results   = args.glm_results
    outpath       = args.outpath
    prefix        = args.prefix
    tdr_key       = args.tdr_key
    num_genes     = args.num_genes
    n_slices      = args.n_slices

    # ------------------------------------------------------------------ #
    # Setup environment                                                   #
    # ------------------------------------------------------------------ #
    os.makedirs(outpath, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Print run information                                               #
    # ------------------------------------------------------------------ #
    print("\n====================== Spateo Interpolation Analysis ========================\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Working directory:", os.getcwd())
    print("Aligned data    :", aligned_data)
    print("PC model        :", pc_model)
    print("Mesh model      :", mesh_model)
    print("Voxel model     :", voxel_model)
    print("GLM results     :", glm_results)
    print("Output directory:", outpath)
    print("Output prefix   :", prefix if prefix else "(none)")
    print("Spateo version  :", st.__version__)
    print("\nKey parameters:")
    print("  tdr_key       :", tdr_key)
    print("  num_genes     :", num_genes)
    print("  n_slices      :", n_slices)
    print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
    print("\n============================================================================\n")

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Step 1: Load data                                                   #
    # ------------------------------------------------------------------ #
    print(">>> Step 1: Loading data\n")
    adata = sc.read_h5ad(aligned_data)
    aligned_pc = st.tdr.read_model(filename=pc_model)
    aligned_mesh = st.tdr.read_model(filename=mesh_model)
    aligned_voxel = st.tdr.read_model(filename=voxel_model)
    glm_data = pd.read_csv(glm_results, index_col=0)
    print("Aligned data: {} cells x {} genes".format(adata.n_obs, adata.n_vars))
    print("GLM results: {} DE genes\n".format(len(glm_data)))

    # ------------------------------------------------------------------ #
    # Step 2: Select genes                                                #
    # ------------------------------------------------------------------ #
    print(">>> Step 2: Selecting top {} genes\n".format(num_genes))
    genes = glm_data.index.tolist()[:num_genes]
    print("Selected genes: {}\n".format(genes))

    # ------------------------------------------------------------------ #
    # Step 3: Raw expression in 3D model                                  #
    # ------------------------------------------------------------------ #
    print(">>> Step 3: Adding raw expression to point cloud model\n")

    pc_index = aligned_pc.point_data["obs_index"].tolist()
    sub = adata[pc_index, genes].X
    if sp.issparse(sub):
        sub = sub.toarray()
    for i, gene_name in enumerate(genes):
        exp = sub[:, i]
        st.tdr.add_model_labels(
            model=aligned_pc,
            labels=exp,
            key_added=gene_name,
            where="point_data",
            inplace=True
        )
    print("Gene expression added to point cloud model\n")

    # ------------------------------------------------------------------ #
    # Step 4: Visualize raw expression                                    #
    # ------------------------------------------------------------------ #
    print(">>> Step 4: Visualizing raw expression\n")

    st.pl.three_d_multi_plot(
        model=aligned_pc,
        key=genes,
        colormap="hot_r",
        opacity=0.5,
        model_style="points",
        jupyter=False,
        text=genes,
        off_screen=True,
        window_size=(800, 800),
        filename=os.path.join(outpath, prefix + "aligned_raw_expr.png"),
        plotter_filename=os.path.join(outpath, prefix + "aligned_raw_expr.html")
    )

    st.pl.three_d_multi_plot(
        model=aligned_pc,
        key=genes,
        colormap="hot_r",
        opacity=0.5,
        model_style="points",
        jupyter=False,
        text=genes,
        off_screen=True,
        window_size=(800, 800),
        filename=os.path.join(outpath, prefix + "aligned_raw_expr.pdf")
    )
    print("Saved: {}aligned_raw_expr.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 5: Gaussian Process interpolation                              #
    # ------------------------------------------------------------------ #
    print(">>> Step 5: Gaussian Process interpolation")
    print("WARNING: This step may take a very long time (minutes to hours)")
    print("         depending on data dimensions and computing resources.\n")

    interpolated_gp_adata = st.tdr.gp_interpolation(
        source_adata=adata.copy(),
        spatial_key=tdr_key,
        keys=genes,
        target_points=np.asarray(aligned_voxel.points)
    )
    print("GP interpolation completed\n")

    # ------------------------------------------------------------------ #
    # Step 6: Construct interpolated point cloud                          #
    # ------------------------------------------------------------------ #
    print(">>> Step 6: Constructing interpolated point cloud\n")

    interpolated_gp_pc, _ = st.tdr.construct_pc(
        adata=interpolated_gp_adata.copy(),
        spatial_key=tdr_key,
        groupby=genes[0]
    )
    _pc_index = interpolated_gp_pc.point_data["obs_index"].tolist()
    for gene_name in genes[1:]:
        _exp = interpolated_gp_adata[_pc_index, gene_name].X.flatten()
        st.tdr.add_model_labels(
            model=interpolated_gp_pc,
            labels=_exp,
            key_added=gene_name,
            where="point_data",
            inplace=True
        )
    print("Interpolated point cloud constructed\n")

    # ------------------------------------------------------------------ #
    # Step 7: Visualize interpolated expression                           #
    # ------------------------------------------------------------------ #
    print(">>> Step 7: Visualizing interpolated expression\n")

    st.pl.three_d_multi_plot(
        model=interpolated_gp_pc,
        key=genes,
        colormap="hot_r",
        opacity=0.5,
        model_style="points",
        jupyter=False,
        off_screen=True,
        text=genes,
        window_size=(800, 800),
        filename=os.path.join(outpath, prefix + "aligned_GP_interpolation.png"),
        plotter_filename=os.path.join(outpath, prefix + "aligned_GP_interpolation.html")
    )

    st.pl.three_d_multi_plot(
        model=interpolated_gp_pc,
        key=genes,
        colormap="hot_r",
        opacity=0.5,
        model_style="points",
        jupyter=False,
        off_screen=True,
        text=genes,
        window_size=(800, 800),
        filename=os.path.join(outpath, prefix + "aligned_GP_interpolation.pdf")
    )
    print("Saved: {}aligned_GP_interpolation.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 8: 3D slice visualization                                      #
    # ------------------------------------------------------------------ #
    print(">>> Step 8: 3D slice visualization\n")

    for gene_name in genes:
        aligned_voxel.point_data[gene_name] = np.asarray(interpolated_gp_adata[:, gene_name].X)

    voxel_slices_x = st.tdr.three_d_slice(
        model=aligned_voxel,
        method="axis",
        n_slices=n_slices,
        axis="x"
    )
    print("Created {} slices along x-axis\n".format(n_slices))

    st.pl.three_d_multi_plot(
        model=st.tdr.collect_models([voxel_slices_x]),
        key=genes,
        model_style="surface",
        colormap="hot_r",
        ambient=0.5,
        jupyter=False,
        off_screen=True,
        shape=(1, 3),
        text=genes,
        window_size=(800, 800),
        filename=os.path.join(outpath, prefix + "aligned_GP_interpolation_slices.png"),
        plotter_filename=os.path.join(outpath, prefix + "aligned_GP_interpolation_slices.html")
    )

    st.pl.three_d_multi_plot(
        model=st.tdr.collect_models([voxel_slices_x]),
        key=genes,
        model_style="surface",
        colormap="hot_r",
        ambient=0.5,
        jupyter=False,
        off_screen=True,
        shape=(1, 3),
        text=genes,
        window_size=(800, 800),
        filename=os.path.join(outpath, prefix + "aligned_GP_interpolation_slices.pdf")
    )
    print("Saved: {}aligned_GP_interpolation_slices.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 9: Save results                                                #
    # ------------------------------------------------------------------ #
    print(">>> Step 9: Saving results\n")

    st.tdr.save_model(
        model=interpolated_gp_pc,
        filename=os.path.join(outpath, prefix + "interpolated_gp_pc.vtk")
    )
    interpolated_gp_adata.write_h5ad(
        os.path.join(outpath, prefix + "interpolated_gp_adata.h5ad"),
        compression="gzip"
    )
    print("Saved: {}interpolated_gp_pc.vtk".format(prefix))
    print("Saved: {}interpolated_gp_adata.h5ad\n".format(prefix))

    elapsed = time.time() - pipeline_start
    print("Total time   : {:.1f} s\n".format(elapsed))

    print("===========================================================================")
    print("Interpolation analysis completed successfully!")
    print("===========================================================================\n")


if __name__ == "__main__":
    main()