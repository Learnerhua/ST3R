#!/path/to/envs/spateo_env/bin/python
import sys, argparse, warnings, os, time, json
from datetime import datetime

# Suppress pkg_resources deprecation warning from spateo
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=UserWarning, module="spateo")

import spateo as st
import pyvista as pv

# Start virtual display for off-screen rendering on servers
pv.start_xvfb()

warnings.filterwarnings("ignore")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate morphological features of 3D reconstructed tissue models."
    )

    # Input/Output arguments
    parser.add_argument("-PC", "--pc_model", required=True, metavar="VTK_FILE",
                        help="Path to point cloud model (vtk file)")
    parser.add_argument("-MS", "--mesh_model", required=True, metavar="VTK_FILE",
                        help="Path to mesh model (vtk file)")
    parser.add_argument("-O", "--outpath", type=str, default=".", metavar="OUTPUT_DIR",
                        help="Output directory (default: current directory)")
    parser.add_argument("-P", "--prefix", type=str, default="", metavar="PREFIX",
                        help="Prefix for output files (default: none)")

    # Key arguments
    parser.add_argument("-KK", "--kde_key", type=str, default="cells_kde", metavar="STR",
                        help="Key for KDE values (default: cells_kde)")

    # Scale factor for morphological calculations
    parser.add_argument("-SF", "--scale_factor", type=float, default=0.5, metavar="FLOAT",
                        help="Scale factor for morphological calculations (default: 0.5)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    pc_model      = args.pc_model
    mesh_model    = args.mesh_model
    outpath       = args.outpath
    prefix        = args.prefix
    kde_key       = args.kde_key
    scale_factor  = args.scale_factor

    # ------------------------------------------------------------------ #
    # Setup environment                                                    #
    # ------------------------------------------------------------------ #
    os.makedirs(outpath, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Print run information                                                #
    # ------------------------------------------------------------------ #
    print("\n====================== Spateo Morphology Analysis ======================\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Working directory:", os.getcwd())
    print("PC model        :", pc_model)
    print("Mesh model      :", mesh_model)
    print("Output directory:", outpath)
    print("Output prefix   :", prefix if prefix else "(none)")
    print("Spateo version  :", st.__version__)
    print("\nKey parameters:")
    print("  kde_key       :", kde_key)
    print("  scale_factor  :", scale_factor)
    print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
    print("\n===================================================================\n")

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Step 1: Load data                                                    #
    # ------------------------------------------------------------------ #
    print(">>> Step 1: Loading models\n")
    aligned_pc = st.tdr.read_model(pc_model)
    aligned_mesh = st.tdr.read_model(mesh_model)
    print("Point cloud model loaded")
    print("Mesh model loaded\n")

    # ------------------------------------------------------------------ #
    # Step 2: Calculate KDE                                                #
    # ------------------------------------------------------------------ #
    print(">>> Step 2: Calculating KDE\n")
    st.tdr.pc_KDE(
        pc=aligned_pc,
        bandwidth=5,
        key_added=kde_key,
        colormap="hot_r",
        inplace=True
    )
    print("KDE calculated\n")

    # ------------------------------------------------------------------ #
    # Step 3: KDE visualization                                            #
    # ------------------------------------------------------------------ #
    print(">>> Step 3: KDE visualization\n")

    st.pl.three_d_plot(
        model=aligned_pc,
        key=kde_key,
        colormap="hot_r",
        opacity=1,
        window_size=(800, 800),
        model_style="points",
        jupyter=False,
        show_axes=True,
        off_screen=True,
        filename=os.path.join(outpath, f"{prefix}aligned_pc_kde.png"),
        plotter_filename=os.path.join(outpath, f"{prefix}aligned_pc_kde.html"),
    )

    st.pl.three_d_plot(
        model=aligned_pc,
        key=kde_key,
        colormap="hot_r",
        opacity=1,
        window_size=(800, 800),
        model_style="points",
        jupyter=False,
        show_axes=True,
        off_screen=True,
        filename=os.path.join(outpath, f"{prefix}aligned_pc_kde.pdf"),
    )
    print("Saved: {}aligned_pc_kde.png(pdf)\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 4: Calculate morphological features                              #
    # ------------------------------------------------------------------ #
    print(">>> Step 4: Calculating morphological features (unit: um)\n")

    # Scale models for morphological calculations
    um_pc_model = aligned_pc.copy()
    um_mesh_model = aligned_mesh.copy()
    um_pc_model.points = um_pc_model.points * scale_factor
    um_mesh_model.points = um_mesh_model.points * scale_factor

    morph = st.tdr.model_morphology(model=um_mesh_model, pc=um_pc_model)

    print("Morphological features:")
    print("  Length (x)    : {:.2f}".format(morph['Length(x)']))
    print("  Width (y)     : {:.2f}".format(morph['Width(y)']))
    print("  Height (z)    : {:.2f}".format(morph['Height(z)']))
    print("  Surface area  : {:.2f}".format(morph['Surface_area']))
    print("  Volume        : {:.2f}".format(morph['Volume']))
    print("  V/SA ratio    : {:.4f}".format(morph['V/SA_ratio']))
    print("  Cell density  : {:.6f}\n".format(morph['cell_density']))

    # ------------------------------------------------------------------ #
    # Step 5: Save results                                                 #
    # ------------------------------------------------------------------ #
    print(">>> Step 5: Saving results\n")

    # Save morphological features as JSON
    with open(os.path.join(outpath, f"{prefix}morph.json"), 'w') as f:
        json.dump(morph, f, indent=4)
    print("Saved: {}morph.json".format(prefix))

    # Save KDE model
    st.tdr.save_model(aligned_pc, os.path.join(outpath, f"{prefix}aligned_pc_KDE_model.vtk"))
    print("Saved: {}aligned_pc_KDE_model.vtk\n".format(prefix))

    elapsed = time.time() - pipeline_start
    print("Total time   : {:.1f} s\n".format(elapsed))

    print("===================================================================")
    print("Morphology analysis completed successfully!")
    print("===================================================================\n")


if __name__ == "__main__":
    main()