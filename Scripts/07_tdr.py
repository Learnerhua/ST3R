#!/oldhome/ouyjh/miniforge3/envs/spateo_env/bin/python
import sys, argparse, warnings, os, time, math
from datetime import datetime

# Suppress pkg_resources deprecation warning from spateo
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=UserWarning, module="spateo")

import torch
import spateo as st
import scanpy as sc
import numpy as np
import pyvista as pv

# Start virtual display for off-screen rendering on servers
pv.start_xvfb()

warnings.filterwarnings("ignore")


def safe_sort_clusters(clusters):
    """Sort clusters, handling both numeric and string cluster names."""
    try:
        return sorted(clusters, key=lambda x: int(x))
    except ValueError:
        return sorted(clusters)


def main():
    parser = argparse.ArgumentParser(
        description="3D reconstruction and visualization of aligned spatial transcriptomics data."
    )

    # Input/Output arguments
    parser.add_argument("-AD", "--aligned_data", required=True, metavar="H5AD_FILE",
                        help="Path to aligned h5ad file")
    parser.add_argument("-RD", "--raw_data", required=True, metavar="H5AD_FILE",
                        help="Path to raw h5ad file")
    parser.add_argument("-O", "--outpath", type=str, default=".", metavar="OUTPUT_DIR",
                        help="Output directory (default: current directory)")
    parser.add_argument("-P", "--prefix", type=str, default="", metavar="PREFIX",
                        help="Prefix for output files (default: none)")

    # Key arguments
    parser.add_argument("-CK", "--cluster_key", type=str, default="squidpy_domains", metavar="STR",
                        help="Key in adata.obs for cluster annotation (default: squidpy_domains)")
    parser.add_argument("-AK", "--aligned_key", type=str, default="aligned_spatial", metavar="STR",
                        help="Key in adata.obsm for aligned spatial coordinates (default: aligned_spatial)")
    parser.add_argument("-TK", "--tdr_key", type=str, default="aligned_spatial_3D", metavar="STR",
                        help="Key to add 3D spatial coordinates (default: aligned_spatial_3D)")
    parser.add_argument("-MK", "--model_key", type=str, default="tissue", metavar="STR",
                        help="Key for the reconstructed model (default: tissue)")

    # Slice distance
    parser.add_argument("-SD", "--slice_distance", type=int, default=30, metavar="INT",
                        help="Slice distance in micrometers (default: 30)")

    # Mesh alpha
    parser.add_argument("-MA", "--mesh_alpha", type=float, default=0.6, metavar="FLOAT",
                        help="Alpha parameter for mesh construction (default: 0.6)")

    # Number of columns for cluster display
    parser.add_argument("-NC", "--ncols", type=int, default=4, metavar="INT",
                        help="Number of columns for cluster display (default: 4)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    aligned_data   = args.aligned_data
    raw_data       = args.raw_data
    outpath        = args.outpath
    prefix         = args.prefix
    cluster_key    = args.cluster_key
    aligned_key    = args.aligned_key
    tdr_key        = args.tdr_key
    model_key      = args.model_key
    slice_distance = args.slice_distance
    mesh_alpha     = args.mesh_alpha
    ncols          = args.ncols

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
    print("\n====================== Spateo 3D Reconstruction ======================\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Working directory:", os.getcwd())
    print("Aligned data    :", aligned_data)
    print("Raw data        :", raw_data)
    print("Output directory:", outpath)
    print("Output prefix   :", prefix if prefix else "(none)")
    print("Device          :", device)
    print("Spateo version  :", st.__version__)
    print("\nKey parameters:")
    print("  cluster_key   :", cluster_key)
    print("  aligned_key   :", aligned_key)
    print("  tdr_key       :", tdr_key)
    print("  model_key     :", model_key)
    print("  slice_distance:", slice_distance, "um")
    print("  mesh_alpha    :", mesh_alpha)
    print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
    print("\n===================================================================\n")

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Step 1: Load data                                                    #
    # ------------------------------------------------------------------ #
    print(">>> Step 1: Loading data\n")
    adata_aligned = sc.read_h5ad(aligned_data)
    adata_raw = sc.read_h5ad(raw_data)
    print("Aligned data: {} cells x {} genes".format(adata_aligned.n_obs, adata_aligned.n_vars))
    print("Raw data    : {} cells x {} genes\n".format(adata_raw.n_obs, adata_raw.n_vars))

    # ------------------------------------------------------------------ #
    # Step 2: Get palette                                                  #
    # ------------------------------------------------------------------ #
    print(">>> Step 2: Getting color palette\n")
    labels = adata_raw.obs[cluster_key].cat.categories.tolist()
    colors = adata_raw.uns[f"{cluster_key}_colors"]
    palette = dict(zip(labels, colors))
    print("Found {} clusters\n".format(len(labels)))

    # ------------------------------------------------------------------ #
    # Step 3: Control slice spacing                                        #
    # ------------------------------------------------------------------ #
    print(">>> Step 3: Setting slice spacing\n")
    xy = adata_aligned.obsm[aligned_key]

    z = (
        adata_aligned.obs['Z']
        .astype('category')
        .cat.codes
        .to_numpy()
        .astype(np.int32)
    )

    slice_spacing = slice_distance * 2
    z_scaled = z * slice_spacing

    adata_aligned.obsm[tdr_key] = np.column_stack([
        xy,
        z_scaled
    ])
    print("Slice distance: {} um\n".format(slice_distance))

    # ------------------------------------------------------------------ #
    # Step 4: Construct point cloud model                                   #
    # ------------------------------------------------------------------ #
    print(">>> Step 4: Constructing point cloud model\n")
    aligned_pc, _ = st.tdr.construct_pc(
        adata=adata_aligned,
        spatial_key=tdr_key,
        groupby=cluster_key,
        key_added=model_key,
        colormap=palette
    )
    print("Point cloud model constructed\n")

    # ------------------------------------------------------------------ #
    # Step 5: 3D visualization                                              #
    # ------------------------------------------------------------------ #
    print(">>> Step 5: 3D visualization\n")

    st.pl.three_d_plot(
        model=aligned_pc,
        key=model_key,
        model_style='points',
        model_size=6,
        show_axes=True,
        jupyter=False,
        off_screen=True,
        window_size=(1200, 1200),
        show_outline=True,
        outline_kwargs={'show_labels': False, 'outline_width': 3},
        filename=os.path.join(outpath, f"{prefix}aligned_pc_3D.png"),
        plotter_filename=os.path.join(outpath, f"{prefix}aligned_pc_3D.html"),
    )

    st.pl.three_d_plot(
        model=aligned_pc,
        key=model_key,
        model_style='points',
        model_size=6,
        show_axes=True,
        jupyter=False,
        off_screen=True,
        window_size=(1200, 1200),
        show_outline=True,
        outline_kwargs={'show_labels': False, 'outline_width': 3},
        filename=os.path.join(outpath, f"{prefix}aligned_pc_3D.pdf"),
    )
    print("Saved: {}aligned_pc_3D.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 6: Three orthogonal views                                        #
    # ------------------------------------------------------------------ #
    print(">>> Step 6: Three orthogonal views\n")

    st.pl.three_d_multi_plot(
        model=st.tdr.collect_models([aligned_pc, aligned_pc, aligned_pc]),
        key=model_key,
        model_style='points',
        model_size=6,
        cpo=['xy', 'xz', 'yz'],
        jupyter=False,
        off_screen=True,
        window_size=(800, 800),
        text=['Spateo xy', 'Spateo xz', 'Spateo yz'],
        filename=os.path.join(outpath, f"{prefix}aligned_pc_3D_multi.png"),
    )

    st.pl.three_d_multi_plot(
        model=st.tdr.collect_models([aligned_pc, aligned_pc, aligned_pc]),
        key=model_key,
        model_style='points',
        model_size=6,
        cpo=['xy', 'xz', 'yz'],
        jupyter=False,
        off_screen=True,
        window_size=(800, 800),
        text=['Spateo xy', 'Spateo xz', 'Spateo yz'],
        filename=os.path.join(outpath, f"{prefix}aligned_pc_3D_multi.pdf"),
    )
    print("Saved: {}aligned_pc_3D_multi.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 7: Clusters display in 3D                                        #
    # ------------------------------------------------------------------ #
    print(">>> Step 7: Clusters display in 3D\n")

    clusters = adata_aligned.obs[cluster_key].astype(str).unique().tolist()
    clusters = safe_sort_clusters(clusters)

    pc_highlight_list = []

    for cluster in clusters:
        highlight_tissues = [cluster]
        pc_highlight = aligned_pc.copy()
        rgba = pc_highlight[model_key + '_rgba']

        # Highlight mask
        mask = adata_aligned.obs[cluster_key].astype(str).isin(highlight_tissues)
        bg_mask = ~mask

        # Non-highlight spots to gray
        rgba[bg_mask, 0] = 180 / 255
        rgba[bg_mask, 1] = 180 / 255
        rgba[bg_mask, 2] = 180 / 255

        # Semi-transparent
        rgba[bg_mask, 3] = 0.05

        # Highlighted tissue keeps original color
        rgba[mask, 3] = 1

        pc_highlight_list.append(pc_highlight)

    nrows = math.ceil(len(clusters) / ncols)
    print("Created {} cluster highlight models (ncols={})\n".format(len(clusters), ncols))

    st.pl.three_d_multi_plot(
        model=st.tdr.collect_models(pc_highlight_list),
        key=model_key,
        model_style='points',
        model_size=4,
        filename=os.path.join(outpath, f"{prefix}aligned_pc_3D_multi_clusters.png"),
        text=["Cluster " + str(i) for i in clusters],
        shape=(nrows, ncols),
        off_screen=True,
        jupyter=False,
    )

    st.pl.three_d_multi_plot(
        model=st.tdr.collect_models(pc_highlight_list),
        key=model_key,
        model_style='points',
        model_size=4,
        filename=os.path.join(outpath, f"{prefix}aligned_pc_3D_multi_clusters.pdf"),
        text=["Cluster " + str(i) for i in clusters],
        shape=(nrows, ncols),
        off_screen=True,
        jupyter=False,
    )
    print("Saved: {}aligned_pc_3D_multi_clusters.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 8: Mesh model                                                    #
    # ------------------------------------------------------------------ #
    print(">>> Step 8: Constructing mesh model\n")

    aligned_mesh, _, _ = st.tdr.construct_surface(
        pc=aligned_pc,
        key_added=model_key,
        alpha=mesh_alpha,
        cs_method="marching_cube",
        cs_args={"mc_scale_factor": 0.8, "dist_sample_num": 20000},
        smooth=2000,
        scale_factor=1.08
    )
    print("Mesh model constructed\n")

    st.pl.three_d_plot(
        model=st.tdr.collect_models([aligned_mesh, aligned_pc]),
        key=model_key,
        model_style=["surface", "points"],
        jupyter=False,
        off_screen=True,
        window_size=(1200, 1200),
        show_outline=False,
        show_axes=True,
        filename=os.path.join(outpath, f"{prefix}aligned_mesh_3D.png"),
        plotter_filename=os.path.join(outpath, f"{prefix}aligned_mesh_3D.html"),
    )

    st.pl.three_d_plot(
        model=st.tdr.collect_models([aligned_mesh, aligned_pc]),
        key=model_key,
        model_style=["surface", "points"],
        jupyter=False,
        off_screen=True,
        window_size=(1200, 1200),
        show_outline=False,
        show_axes=True,
        filename=os.path.join(outpath, f"{prefix}aligned_mesh_3D.pdf"),
    )
    print("Saved: {}aligned_mesh_3D.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 9: Voxel model                                                   #
    # ------------------------------------------------------------------ #
    print(">>> Step 9: Constructing voxel model\n")

    aligned_voxel, _ = st.tdr.voxelize_mesh(
        mesh=aligned_mesh,
        voxel_pc=None,
        key_added=model_key,
        label="voxel",
        color="gainsboro",
        smooth=500
    )
    print("Voxel model constructed\n")

    st.pl.three_d_plot(
        model=aligned_voxel,
        key=model_key,
        off_screen=True,
        window_size=(1200, 1200),
        show_outline=False,
        show_axes=True,
        filename=os.path.join(outpath, f"{prefix}aligned_voxel_3D.png"),
        plotter_filename=os.path.join(outpath, f"{prefix}aligned_voxel_3D.html"),
    )

    st.pl.three_d_plot(
        model=aligned_voxel,
        key=model_key,
        off_screen=True,
        window_size=(1200, 1200),
        show_outline=False,
        show_axes=True,
        filename=os.path.join(outpath, f"{prefix}aligned_voxel_3D.pdf"),
    )
    print("Saved: {}aligned_voxel_3D.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 10: Save all models                                              #
    # ------------------------------------------------------------------ #
    print(">>> Step 10: Saving all models\n")

    st.tdr.save_model(
        model=aligned_pc,
        filename=os.path.join(outpath, f"{prefix}aligned_pc_model.vtk")
    )
    st.tdr.save_model(
        model=aligned_mesh,
        filename=os.path.join(outpath, f"{prefix}aligned_mesh_model.vtk")
    )
    st.tdr.save_model(
        model=aligned_voxel,
        filename=os.path.join(outpath, f"{prefix}aligned_voxel_model.vtk")
    )
    adata_aligned.write_h5ad(os.path.join(outpath, f"{prefix}tdr.h5ad"), compression='gzip')
    print("Saved: {}aligned_pc_model.vtk".format(prefix))
    print("Saved: {}aligned_mesh_model.vtk".format(prefix))
    print("Saved: {}aligned_voxel_model.vtk".format(prefix))
    print("Saved: {}tdr.h5ad\n".format(prefix))

    elapsed = time.time() - pipeline_start
    print("Total time   : {:.1f} s\n".format(elapsed))

    print("===================================================================")
    print("3D reconstruction completed successfully!")
    print("===================================================================\n")


if __name__ == "__main__":
    main()