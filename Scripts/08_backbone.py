#!/oldhome/ouyjh/miniforge3/envs/spateo_env/bin/python
import sys, argparse, warnings, os, time
from datetime import datetime

# Suppress pkg_resources deprecation warning from spateo
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=UserWarning, module="spateo")

import numpy as np
import spateo as st
import scanpy as sc
import pyvista as pv

# Start virtual display for off-screen rendering on servers
pv.start_xvfb()

warnings.filterwarnings("ignore")


def main():
    parser = argparse.ArgumentParser(
        description="Backbone construction and differential expression analysis for 3D spatial transcriptomics data."
    )

    # Input/Output arguments
    parser.add_argument("-AD", "--aligned_data", required=True, metavar="H5AD_FILE",
                        help="Path to aligned h5ad file from TDR step")
    parser.add_argument("-PC", "--pc_model", required=True, metavar="VTK_FILE",
                        help="Path to point cloud model VTK file")
    parser.add_argument("-MS", "--mesh_model", required=True, metavar="VTK_FILE",
                        help="Path to mesh model VTK file")
    parser.add_argument("-O", "--outpath", type=str, default=".", metavar="OUTPUT_DIR",
                        help="Output directory (default: current directory)")
    parser.add_argument("-P", "--prefix", type=str, default="", metavar="PREFIX",
                        help="Prefix for output files (default: none)")

    # Backbone parameters
    parser.add_argument("-NN", "--num_nodes", type=int, default=30, metavar="INT",
                        help="Number of backbone nodes (default: 30)")

    # DE analysis parameters
    parser.add_argument("-DK", "--de_key", type=str, default="glm_degs", metavar="STR",
                        help="Key for storing DE results in adata.uns (default: glm_degs)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    aligned_data = args.aligned_data
    pc_model     = args.pc_model
    mesh_model   = args.mesh_model
    outpath      = args.outpath
    prefix       = args.prefix
    num_nodes    = args.num_nodes
    de_key       = args.de_key

    # ------------------------------------------------------------------ #
    # Setup environment                                                   #
    # ------------------------------------------------------------------ #
    os.makedirs(outpath, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Print run information                                               #
    # ------------------------------------------------------------------ #
    print("\n====================== Spateo Backbone Analysis ======================\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Working directory:", os.getcwd())
    print("Aligned data    :", aligned_data)
    print("PC model        :", pc_model)
    print("Mesh model      :", mesh_model)
    print("Output directory:", outpath)
    print("Output prefix   :", prefix if prefix else "(none)")
    print("Spateo version  :", st.__version__)
    print("\nKey parameters:")
    print("  num_nodes     :", num_nodes)
    print("  de_key        :", de_key)
    print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
    print("\n===================================================================\n")

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Step 1: Load data                                                   #
    # ------------------------------------------------------------------ #
    print(">>> Step 1: Loading data\n")
    adata = sc.read_h5ad(aligned_data)
    aligned_pc = st.tdr.read_model(filename=pc_model)
    aligned_mesh = st.tdr.read_model(filename=mesh_model)
    print("Aligned data: {} cells x {} genes\n".format(adata.n_obs, adata.n_vars))

    # ------------------------------------------------------------------ #
    # Step 2: Construct backbone                                          #
    # ------------------------------------------------------------------ #
    print(">>> Step 2: Constructing backbone\n")

    backbone, backbone_length, _ = st.tdr.construct_backbone(
        model=aligned_pc,
        rd_method="ElPiGraph",
        num_nodes=num_nodes,
        color="gainsboro"
    )
    print("Backbone length: {:.2f}".format(backbone_length))

    # Check backbone coordinates
    pp = np.asarray(backbone.points)
    for i, axis in enumerate(['X', 'Y', 'Z']):
        print("  {}: {:.2f} ~ {:.2f}  range: {:.2f}".format(
            axis, pp[:, i].min(), pp[:, i].max(), pp[:, i].max() - pp[:, i].min()
        ))
    print()

    # ------------------------------------------------------------------ #
    # Step 3: Visualize backbone with models                              #
    # ------------------------------------------------------------------ #
    print(">>> Step 3: Visualizing backbone with models\n")

    st.pl.three_d_plot(
        model=st.tdr.collect_models([aligned_pc, aligned_mesh, backbone]),
        key="backbone",
        colormap=["gainsboro", "gainsboro", "orangered"],
        model_style=["points", "surface", "wireframe"],
        model_size=[3, None, 8],
        opacity=[0.3, 0.3, 1.0],
        show_legend=False,
        jupyter=False,
        off_screen=True,
        window_size=(800, 800),
        show_outline=False,
        show_axes=True,
        filename=os.path.join(outpath, prefix + 'backbone_3D.png'),
        plotter_filename=os.path.join(outpath, prefix + 'backbone_3D.html')
    )

    st.pl.three_d_plot(
        model=st.tdr.collect_models([aligned_pc, aligned_mesh, backbone]),
        key="backbone",
        colormap=["gainsboro", "gainsboro", "orangered"],
        model_style=["points", "surface", "wireframe"],
        model_size=[3, None, 8],
        opacity=[0.3, 0.3, 1.0],
        show_legend=False,
        jupyter=False,
        off_screen=True,
        window_size=(800, 800),
        show_outline=False,
        show_axes=True,
        filename=os.path.join(outpath, prefix + 'backbone_3D.pdf')
    )
    print("Saved: {}backbone_3D.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 4: Map points to backbone                                      #
    # ------------------------------------------------------------------ #
    print(">>> Step 4: Mapping points to backbone\n")

    st.tdr.map_points_to_backbone(
        model=aligned_pc,
        backbone_model=backbone,
        nodes_key="nodes",
        key_added="backbone",
        inplace=True,
    )

    print("Total number of cells in adata:", adata.n_obs)
    print("Number of backbone nodes:", len(np.unique(aligned_pc.point_data["backbone"])))
    adata.obs["backbone"] = aligned_pc.point_data["backbone"]
    print("Mapping to adata successful, number of cells:", adata.n_obs)
    print("Cell count per node:\n", adata.obs["backbone"].value_counts().sort_index())
    print()

    # ------------------------------------------------------------------ #
    # Step 5: Visualize backbone area                                     #
    # ------------------------------------------------------------------ #
    print(">>> Step 5: Visualizing backbone area\n")

    st.pl.three_d_plot(
        model=st.tdr.collect_models([aligned_pc, backbone]),
        key="backbone",
        colormap="rainbow",
        model_style=["points", "wireframe"],
        model_size=[2, 8],
        opacity=[0.5, 0.8],
        show_legend=False,
        jupyter=False,
        off_screen=True,
        window_size=(800, 800),
        show_outline=False,
        show_axes=True,
        filename=os.path.join(outpath, prefix + 'backbone_area.png'),
        plotter_filename=os.path.join(outpath, prefix + 'backbone_area.html')
    )

    st.pl.three_d_plot(
        model=st.tdr.collect_models([aligned_pc, backbone]),
        key="backbone",
        colormap="rainbow",
        model_style=["points", "wireframe"],
        model_size=[2, 8],
        opacity=[0.5, 0.8],
        show_legend=False,
        jupyter=False,
        off_screen=True,
        window_size=(800, 800),
        show_outline=False,
        show_axes=True,
        filename=os.path.join(outpath, prefix + 'backbone_area.pdf')
    )
    print("Saved: {}backbone_area.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 6: Differential expression analysis                            #
    # ------------------------------------------------------------------ #
    print(">>> Step 6: Differential expression analysis")
    print("WARNING: This step may take a very long time (hours to days)")
    print("         depending on data dimensions and computing resources.\n")

    # Suppress statsmodels warnings
    warnings.filterwarnings("ignore", category=UserWarning, module="statsmodels")

    st.tl.glm_degs(
        adata=adata,
        fullModelFormulaStr='~cr(backbone, df=3)',
        key_added=de_key,
        qval_threshold=0.05,
        llf_threshold=-1000,
    )

    glm_data = adata.uns[de_key]["glm_result"]
    glm_data.to_csv(os.path.join(outpath, "glm_data.csv"))
    print("Saved: glm_data.csv")
    print("Number of DE genes: {}\n".format(len(glm_data)))

    # ------------------------------------------------------------------ #
    # Step 7: Visualize top genes glm fit                                 #
    # ------------------------------------------------------------------ #
    print(">>> Step 7: Visualizing top genes glm fit\n")

    prefix_tmp = prefix.replace("_", "")

    st.pl.glm_fit(
        adata=adata,
        genes=glm_data.index.tolist()[:9],
        ncols=3,
        feature_x="backbone",
        feature_y="expression",
        glm_key=de_key,
        save_show_or_return="save",
        save_kwargs={
            "path": os.path.join(outpath, "top9Genes_glm_fit"),
            "prefix": prefix_tmp,
            "dpi": 300,
            "ext": "png",
            "transparent": False
        }
    )

    st.pl.glm_fit(
        adata=adata,
        genes=glm_data.index.tolist()[:9],
        ncols=3,
        feature_x="backbone",
        feature_y="expression",
        glm_key=de_key,
        save_show_or_return="save",
        save_kwargs={
            "path": os.path.join(outpath, "top9Genes_glm_fit"),
            "prefix": prefix_tmp,
            "ext": "pdf",
            "transparent": False
        }
    )
    print("Saved: {}top9Gene_glm_fit.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 8: Save results                                                #
    # ------------------------------------------------------------------ #
    print(">>> Step 8: Saving results\n")

    st.tdr.save_model(
        model=backbone,
        filename=os.path.join(outpath, prefix + "backbone_model.vtk")
    )
    adata.write_h5ad(os.path.join(outpath, prefix + "backbone.h5ad"), compression='gzip')
    print("Saved: {}backbone_model.vtk".format(prefix))
    print("Saved: {}backbone.h5ad\n".format(prefix))

    elapsed = time.time() - pipeline_start
    print("Total time   : {:.1f} s\n".format(elapsed))

    print("===================================================================")
    print("Backbone analysis completed successfully!")
    print("===================================================================\n")


if __name__ == "__main__":
    main()