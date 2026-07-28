# ST3R

**ST3R** — **S**tereo-seq **3D** **T**issue **R**econstruction Pipeline

A focused 3D reconstruction pipeline for Stereo-seq spatial transcriptomics data, performing multi-slice alignment, 3D tissue modeling, backbone extraction, morphological quantification, and Gaussian Process gene expression interpolation on a dense 3D voxel grid, powered by [Spateo](https://github.com/aristoteleo/spateo-release).

> **中文文档**: [README.zh-CN.md](README.zh-CN.md) · [安装指南中文版](INSTALL.zh-CN.md)

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](#license)
[![Python 3.8 / 3.10 / 3.12](https://img.shields.io/badge/Python-3.8%20%2F%203.10%20%2F%203.12-blue.svg)](https://www.python.org)
[![Platform: Linux](https://img.shields.io/badge/Platform-Linux-lightgrey.svg)](#)
[![Powered by Spateo](https://img.shields.io/badge/Powered%20by-Spateo-00A1DE.svg)](https://github.com/aristoteleo/spateo-release)

[![Stereo-seq](https://img.shields.io/badge/Compatible-Stereo--seq-FF6B6B)](#)
[![GPU Accelerated](https://img.shields.io/badge/GPU-CUDA%20Optional-76B900.svg)](#installation)
[![Steps: 11](https://img.shields.io/badge/Pipeline-11%20Steps-orange.svg)](#pipeline-architecture)

</div>

---

## Overview

**ST3R** is a 3D-reconstruction-focused pipeline for Stereo-seq spatial transcriptomics data. It centers on Spateo's `tdr` (Three-Dimensional Reconstruction) module to transform multi-slice 2D Stereo-seq data into 3D point cloud, mesh, and voxel tissue models, then quantifies morphology and interpolates gene expression into continuous 3D space. The pipeline ships a 4-step upstream preprocessing chain (Steps 01-06) to produce TDR-ready input and a downstream reporting layer (Step 11) that packages every artifact into a client-deliverable HTML/PDF report. Suitable for Stereo-seq tissue slices requiring 3D coordinate reconstruction, developmental biology studies, organoid 3D modeling, and spatially-resolved transcriptomic analysis on reconstructed anatomy.

---

## Key Features

- **3D Reconstruction as the Core**: Centered on Spateo's `tdr` module — Steps 07-10 cover point cloud / mesh / voxel construction, backbone extraction, morphological quantification, and GP-based gene expression interpolation.
- **GPU-Accelerated 3D Pipeline**: Steps 06 (alignment) and 07 (3D reconstruction) support CUDA acceleration, with automatic CPU fallback.
- **Standardized H5AD I/O**: All intermediate artifacts use the AnnData format for seamless Scanpy/Squidpy integration.
- **Built-in Spatial Domain Detection**: Squidpy-based joint graph clustering identifies tissue domains prior to alignment.
- **Multi-Model 3D Output**: Generates point cloud, mesh surface, and voxel models in VTK format for ParaView / PyVista.
- **Backbone Differential Expression**: GLM-based DE analysis along the principal anatomical axis of the reconstructed tissue.
- **Gaussian Process Interpolation**: Predicts gene expression on a dense 3D voxel grid from sparse spatial measurements.
- **Automated Report Generation**: Step 11 aggregates all PNG/HTML/CSV/JSON outputs into interactive HTML (~47 KB) and paginated PDF (~40 MB) reports — no re-analysis required.
- **Flexible Report Variants**: Two report scripts available — `11_report.py` (full pipeline, 7 chapters) and `11_report_subset.py` (alignment + TDR only).

---

## Pipeline Architecture

The pipeline consists of 11 sequential Python scripts organized in three layers around the 3D reconstruction core:

| Layer | Steps | Purpose |
|-------|-------|---------|
| **Upstream preprocessing** | 01–05 | GEF → H5AD → QC → spatial domains → clean input |
| **3D reconstruction core** ⚡ | 06–10 | Alignment → 3D models → backbone → morphology → GP interpolation |
| **Reporting** | 11 | Aggregate artifacts into HTML/PDF reports |

| Step | Script | Conda Env | Input | Output | Key Tools |
|------|--------|-----------|-------|--------|-----------|
| 01 | `01_gef2h5ad.py` | `stereopy` | GEF file(s) / TSV config | `{sample}.h5ad` | stereo |
| 02 | `02_concat.py` | `scanpy` | H5AD directory | `*_concated.h5ad` | anndata |
| 03 | `03_preprocess.py` | `scanpy` | Concatenated H5AD | `*_preprocessed.h5ad` + QC plots | scanpy, harmonypy, scrublet |
| 04 | `04_squidpy.py` | `scanpy` | Preprocessed H5AD | `*_squidpy.h5ad` + spatial domain plots | squidpy, scanpy |
| 05 | `05_dataConvert.py` | `scanpy` | Squidpy H5AD | `*_compatible.h5ad` | scanpy |
| 06 | `06_align.py` ⚡ | `spateo_env` | Compatible H5AD | `*_adata_aligned.h5ad` + alignment plots | spateo, torch (CUDA) |
| **07** | **`07_tdr.py` ⚡** | `spateo_env` | Aligned H5AD + compatible H5AD | `*.vtk` (PC/mesh/voxel) + `*_tdr.h5ad` | **spateo.tdr**, pyvista, torch (CUDA) |
| **08** | **`08_backbone.py`** | `spateo_env` | TDR H5AD + VTK models | `glm_data.csv` + backbone model | **spateo.tdr**, pyvista |
| **09** | **`09_morph.py`** | `spateo_env` | PC + mesh VTK models | `*_morph.json` + KDE plot | **spateo.tdr**, pyvista |
| **10** | **`10_interpolation.py`** | `spateo_env` | TDR H5AD + VTK models + GLM CSV | `*_interpolated_gp_adata.h5ad` + plots | **spateo.tdr.gp_interpolation** |
| 11 | `11_report.py` | `scanpy` | Output root (Steps 01-10) | `report.html` + `report.pdf` + assets | jinja2, playwright, pypdf2 |

> ⚡ = GPU-accelerated (CUDA optional; falls back to CPU automatically)
> **Bold rows = the 3D reconstruction core** (Steps 07-10), all built on Spateo's `tdr` module.
> Step 11 shares the `scanpy` environment (Steps 02-05) — its dependencies (jinja2, playwright, pypdf2) are co-installed.
> Step 11 does not re-run any analysis — it only aggregates the artifacts produced by Steps 01-10.

---

## Requirements

### Required System Tools

- **Python** (v3.10+) — runtime
- **conda** or **mamba** — environment management
- **CUDA Toolkit** (optional, v11.8+) — GPU acceleration for alignment and 3D reconstruction
- **Xvfb** — off-screen rendering for PyVista on headless servers

### Required Python Packages

```bash
pip install stereo spateo scanpy squidpy anndata pyvista
pip install harmonypy scrublet scipy pandas numpy
```

### Input Data Convention

- **GEF file format**: Stereo-seq cell-segmented (`*.cellbin.gef`) or fixed-bin (`*.gef`) expression files
- **Sample list TSV**: Three required columns — `gef_path` (absolute path to GEF), `sample_id` (unique sample identifier), `slice_id` (slice label used for sorting; numeric suffix drives Z-axis ordering)
- **Bin type**: `cell_bins` recommended (cell-segmented); `bins` for fixed-resolution data
- **Minimum slices**: At least 2 spatial slices required for downstream alignment

---

## Installation

> ⚠️ **Note**: Installation commands below are **reference examples only**. For complete step-by-step instructions, see **[INSTALL.md](INSTALL.md)**.

The pipeline uses 3 conda environments to satisfy strict Python-version constraints:

| Conda env | Python | Steps | Key packages |
|-----------|--------|-------|--------------|
| `stereopy` | 3.8 | 01 | stereo, pandas |
| `scanpy` | 3.12 | 02-05, 11 | scanpy, squidpy, anndata, jinja2, playwright, pypdf2 |
| `spateo_env` | 3.10 | 06-10 | spateo, torch (CUDA), pyvista |

### Quick Install

```bash
mamba env create -f envs/stereopy.yml
mamba env create -f envs/scanpy.yml
mamba env create -f envs/spateo_env.yml

# Step 11 requires headless Chromium for PDF rendering
conda run -n scanpy playwright install chromium
```

### Next Steps

After installation, use the bundled `run.sh` wrapper (no setup needed — just `chmod +x run.sh`) to invoke scripts without modifying their hardcoded shebangs:

```bash
./run.sh 01_gef2h5ad.py -C rawData/sample_list.tsv -BT cell_bins -O Output/01_gef2h5ad
```

For GPU setup, mirror configuration, troubleshooting, and detailed verification commands, see **[INSTALL.md](INSTALL.md)**.

---

## Usage

The pipeline supports two execution modes. Most users should start with **Mode 1**.

### Mode 1 — Agent-driven Execution (Recommended)

This project is designed to be driven by an AI coding agent (e.g., Claude Code). Each script is self-documenting via its `--help` flag and follows a consistent CLI pattern (`-I/--input`, `-O/--outpath`, `-P/--prefix`).

**Setup:**

1. Open the project root in your AI coding agent.
2. The agent reads the scripts and walks through every parameter.
3. Issue a natural-language request — the agent validates the environment and confirms each step before execution.

**Example prompts:**

```
Run the full Spateo pipeline on a sample list TSV with prefix "Sol_".
```

```
Re-run only Steps 08-10 using the existing VTK models in Sol_test/07_tdr/.
```

### Mode 2 — Manual Step-by-step Execution

> ⚠️ **Note**: The shell commands below are **reference examples only**. Parameter values (especially QC thresholds `-minG`, `-maxG`, `-minU`, `-maxU`, `-maxMT`, `-maxHB`) must be tuned to your specific dataset. Always run `python <script>.py --help` first to view the full parameter list.

#### Quick Start — Run the Full Pipeline

```bash
# Step 01: Convert GEF to H5AD
./Scripts/01_gef2h5ad.py -C rawData/sample_list.tsv -BT cell_bins -O Output/01_gef2h5ad

# Step 02: Concatenate H5AD files
./Scripts/02_concat.py -I Output/01_gef2h5ad -O Output/02_concat -P Sol_

# Step 03: Preprocess (set QC thresholds appropriate for your data)
./Scripts/03_preprocess.py -I Output/02_concat/Sol_concated.h5ad -BK slice_id -P Sol_ -O Output/03_preprocess \
    -minG 200 -maxG 2000 -minU 200 -maxU 6000 -minC 3 -maxMT 5 -maxHB 5

# Step 04: Spatial domain detection
./Scripts/04_squidpy.py -I Output/03_preprocess/Sol_preprocessed.h5ad -LK slice_id -P Sol_ -O Output/04_squidpy \
    -R 1.2 -WS 0.8

# Step 05: Data cleaning for Spateo compatibility
./Scripts/05_dataConvert.py -I Output/04_squidpy/Sol_squidpy.h5ad -P Sol_ -O Output/05_dataConvert

# Step 06: 3D alignment (GPU-accelerated)
./Scripts/06_align.py -I Output/05_dataConvert/Sol_compatible.h5ad -P Sol_ -O Output/06_alignment

# Step 07: 3D tissue reconstruction (GPU-accelerated)
./Scripts/07_tdr.py -AD Output/06_alignment/Sol_adata_aligned.h5ad -RD Output/05_dataConvert/Sol_compatible.h5ad -P Sol_ -O Output/07_tdr

# Step 08: Backbone extraction + GLM DE
./Scripts/08_backbone.py -AD Output/07_tdr/Sol_tdr.h5ad -PC Output/07_tdr/Sol_aligned_pc_model.vtk \
    -MS Output/07_tdr/Sol_aligned_mesh_model.vtk -P Sol_ -O Output/08_backbone

# Step 09: Morphological features
./Scripts/09_morph.py -PC Output/07_tdr/Sol_aligned_pc_model.vtk -MS Output/07_tdr/Sol_aligned_mesh_model.vtk -P Sol_ -O Output/09_Morph

# Step 10: Gaussian process interpolation
./Scripts/10_interpolation.py -AD Output/07_tdr/Sol_tdr.h5ad -PC Output/07_tdr/Sol_aligned_pc_model.vtk \
    -MS Output/07_tdr/Sol_aligned_mesh_model.vtk -VX Output/07_tdr/Sol_aligned_voxel_model.vtk \
    -GL Output/08_backbone/glm_data.csv -P Sol_ -O Output/10_interpolation -NG 3 -NS 15

# Step 11a: Full client report (HTML + PDF, all 7 chapters)
./Scripts/11_report.py -I Output -P Sol_ -SL rawData/sample_list.tsv -O Output/11_report

# Step 11b: Subset report (alignment + TDR only, 2 chapters)
./Scripts/11_report_subset.py -I Output -P Sol_ -SL rawData/sample_list.tsv -O Output/11_report_subset
```

For complete parameter documentation for every script, see [`Spateo_pipeline_SOP.pdf`](Spateo_pipeline_SOP.pdf).

#### Single-step Debugging

Each script can be re-run independently given the required inputs from prior steps. Use the table in **Pipeline Architecture** to identify dependencies.

---

## Output Files

Each step writes its outputs to a numbered subdirectory under `Output/`:

| Step | Output Directory | Key Files |
|------|------------------|-----------|
| 01 | `Output/01_gef2h5ad/` | `{sample_id}.cellbin.h5ad` (one per sample) |
| 02 | `Output/02_concat/` | `{prefix}_concated.h5ad` |
| 03 | `Output/03_preprocess/` | `{prefix}_preprocessed.h5ad`, QC violin/scatter plots, PCA, UMAP, marker gene heatmap, `{prefix}_all_markers.csv` |
| 04 | `Output/04_squidpy/` | `{prefix}_squidpy.h5ad`, spatial domain grids, leiden comparison UMAP |
| 05 | `Output/05_dataConvert/` | `{prefix}_compatible.h5ad` |
| 06 | `Output/06_alignment/` | `{prefix}_adata_aligned.h5ad`, before/after alignment slice plots, overlap comparison |
| 07 | `Output/07_tdr/` | `{prefix}_aligned_pc_model.vtk`, `{prefix}_aligned_mesh_model.vtk`, `{prefix}_aligned_voxel_model.vtk`, `{prefix}_tdr.h5ad`, interactive HTML viewers |
| 08 | `Output/08_backbone/` | `{prefix}_backbone_model.vtk`, `glm_data.csv`, top-9 DE gene GLM fit plots, `{prefix}_backbone.h5ad` |
| 09 | `Output/09_Morph/` | `{prefix}_morph.json` (length/width/height/SA/volume/density), KDE heatmap, `{prefix}_aligned_pc_KDE_model.vtk` |
| 10 | `Output/10_interpolation/` | `{prefix}_interpolated_gp_adata.h5ad`, raw vs interpolated expression plots, 3D slice views, `{prefix}_interpolated_gp_pc.vtk` |
| 11 | `Output/11_report/` | `report.html` (~47 KB interactive), `report.pdf` (~40 MB paginated), `assets/`, `images/` (24 figures in 7 chapter subdirs), `viz/` (interactive 3D iframes), `data/` (CSV/JSON) |

---

## Example: Reference Test Run (33 Slices)

A complete test run on 33 Stereo-seq slices is included in `Sol_test/`. Example outputs:

### 3.1 Step 03: Preprocessing (QC, Normalization, Clustering)

**Before vs After QC**:

| Before QC | After QC |
|-----------|----------|
| ![before](Sol_test/03_preprocess/Sol_before_QC_violin.png) | ![after](Sol_test/03_preprocess/Sol_after_QC_violin.png) |

**PCA Variance & UMAP Clustering**:

| PCA Variance Ratio | UMAP (Leiden) |
|--------------------|---------------|
| ![pca_var](Sol_test/03_preprocess/Sol_pca_variance_ratio.png) | ![umap](Sol_test/03_preprocess/Sol_umap_leiden.png) |

---

### 3.2 Step 04: Spatial Domain Detection (Squidpy)

**Spatial Domain Grid**:

![domains_grid](Sol_test/04_squidpy/Sol_squidpy_domains_grid.png)

**Leiden vs Squidpy Domain Comparison**:

![leiden_squidpy](Sol_test/04_squidpy/Sol_leiden_squidpy.png)

---

### 3.3 Step 06: 3D Alignment

**Before Alignment**:

![before](Sol_test/06_alignment/Sol_squidpy_2Dslices.png)

**After Alignment**:

![after](Sol_test/06_alignment/Sol_aligned_2Dslices.png)

**Before/After Overlap Comparison**:

![overlap](Sol_test/06_alignment/Sol_aligned_2Dslices_overlap.png)

---

### 3.4 Step 07: 3D Tissue Reconstruction (TDR)

**Point Cloud, Mesh, and Voxel Models**:

| Point Cloud | Mesh Surface | Voxel |
|-------------|--------------|-------|
| ![pc](Sol_test/07_tdr/Sol_aligned_pc_3D.png) | ![mesh](Sol_test/07_tdr/Sol_aligned_mesh_3D.png) | ![voxel](Sol_test/07_tdr/Sol_aligned_voxel_3D.png) |

**Orthogonal Three-View Projection**:

![multi](Sol_test/07_tdr/Sol_aligned_pc_3D_multi.png)

---

### 3.5 Step 08: Backbone Extraction + GLM Differential Expression

| Backbone 3D | Backbone Node Coloring |
|-------------|------------------------|
| ![bb3d](Sol_test/08_backbone/Sol_backbone_3D.png) | ![bb_area](Sol_test/08_backbone/Sol_backbone_area.png) |

**Top 9 DE Genes — GLM Fit**:

![glm_fit](Sol_test/08_backbone/Sol_top9Genes_glm_fit.png)

---

### 3.6 Step 09: Morphological Features

**Cell Density KDE Heatmap**:

![kde](Sol_test/09_Morph/Sol_aligned_pc_kde.png)

---

### 3.7 Step 10: Gaussian Process Interpolation

**Raw Expression**:

![raw](Sol_test/10_interpolation/Sol_aligned_raw_expr.png)

**GP-Interpolated Expression**:

![gp](Sol_test/10_interpolation/Sol_aligned_GP_interpolation.png)

**3D Slice Views**:

![slices](Sol_test/10_interpolation/Sol_aligned_GP_interpolation_slices.png)

---

## Notes

1. **Step 08 GLM runtime**: The GLM differential expression step scales with cell count and can take hours to days on large datasets. Consider reducing `-NN` (number of backbone nodes) for exploratory runs.
2. **Step 10 GP interpolation runtime**: Gaussian Process interpolation on dense voxel grids is computationally expensive. The `-NG` (number of genes) and `-NS` (number of slices) parameters directly control memory and runtime.
3. **Step 11 PDF size**: The full Step 11a PDF is ~40 MB because 3D models are embedded as static screenshots. Use `--no-pdf` to generate only the lightweight HTML (~47 KB) for quick preview.
4. **Naming convention coupling**: The `{prefix}` (e.g., `Sol_`) used in Steps 02-10 must remain consistent. If you change it, every downstream `-I` path and `-P` flag must be updated in lockstep.
5. **Sample list TSV paths**: All `gef_path` entries in `sample_list.tsv` must be absolute paths. Relative paths will fail silently during batch processing.
6. **Bin type choice**: `cell_bins` produces cell-resolved expression data and is recommended for most use cases. `bins` mode requires choosing `-BS` (default 50, in microns) and produces fixed-grid data better suited for cell-type deconvolution.
7. **Z-axis ordering**: Slice Z-coordinates in the 3D model are derived from the numeric suffix of `slice_id` (e.g., `Sol_1`, `Sol_2`, ...). Use a consistent numeric suffix scheme.
8. **GPU memory**: For datasets with > 500K cells, ensure at least 16 GB of GPU memory is available for Steps 06 and 07.
9. **Intermediate disk space**: A 33-slice cell-bin dataset produces ~5-10 GB of intermediate H5AD/VTK files. Plan disk accordingly; do not delete intermediate outputs between steps unless you intend to re-run.

---

## Project Structure

```
.
├── README.md                       # This file
├── INSTALL.md                      # Detailed installation guide
├── Spateo_pipeline_SOP.pdf         # Detailed standard operating procedure (PDF)
├── run.sh                          # Conda-aware wrapper to invoke any step
├── Scripts/                        # Pipeline scripts (11 sequential steps)
│   ├── 01_gef2h5ad.py              # Step 01: GEF → H5AD conversion
│   ├── 02_concat.py                # Step 02: H5AD concatenation
│   ├── 03_preprocess.py            # Step 03: QC + normalization + clustering
│   ├── 04_squidpy.py               # Step 04: Spatial domain detection
│   ├── 05_dataConvert.py           # Step 05: H5AD cleaning
│   ├── 06_align.py                 # Step 06: 3D alignment (GPU)
│   ├── 07_tdr.py                   # Step 07: 3D tissue reconstruction (GPU)
│   ├── 08_backbone.py              # Step 08: Backbone + GLM DE
│   ├── 09_morph.py                 # Step 09: Morphological features
│   ├── 10_interpolation.py         # Step 10: GP interpolation
│   ├── 11_report.py                # Step 11a: Full client report (HTML + PDF)
│   ├── 11_report_subset.py         # Step 11b: Subset report (alignment + TDR only)
│   └── Report_config/              # Step 11 configuration package
│       ├── config.yaml             # Report section/chapter toggles
│       ├── templates/              # Jinja2 HTML template + PDF cover
│       ├── lib/                    # Render, content, and collect helpers
│       └── assets/                 # CSS, logos, workflow diagram
├── envs/                           # Conda environment definitions
│   ├── stereopy.yml
│   ├── scanpy.yml
│   └── spateo_env.yml
└── Sol_test/                       # Reference test run outputs (33 slices)
    ├── 03_preprocess/              # QC plots, UMAP, markers
    ├── 04_squidpy/                 # Spatial domain grids
    ├── 06_alignment/               # Before/after alignment plots
    ├── 07_tdr/                     # 3D reconstruction outputs
    ├── 08_backbone/                # Backbone + GLM DE results
    ├── 09_Morph/                   # Morphology outputs
    └── 10_interpolation/           # GP interpolation outputs
```

---

## Contact

For questions, bug reports, or feature requests, please submit an Issue or Pull Request.

**Email**: oyjh417701@163.com

## Copyright

Copyright (c) 2026 OYJH. All Rights Reserved.

## License

MIT License

## Acknowledgments

This pipeline integrates the following excellent open-source tools:

- [Spateo](https://github.com/aristoteleo/spateo-release) — spatial transcriptomics alignment, 3D reconstruction, GP interpolation
- [StereoPy](https://github.com/BGI-Shenzhen/StereoPy) — Stereo-seq GEF file I/O
- [Scanpy](https://github.com/scverse/scanpy) — single-cell preprocessing and analysis
- [Squidpy](https://github.com/scverse/squidpy) — spatial neighborhood graph construction
- [PyVista](https://github.com/pyvista/pyvista) — 3D mesh and volume rendering
- [AnnData](https://github.com/scverse/anndata) — annotated data matrix format
- [Jinja2](https://github.com/pallets/jinja) — report templating
- [Playwright](https://github.com/microsoft/playwright) — headless Chromium HTML-to-PDF rendering
- [PyPDF2](https://github.com/py-pdf/PyPDF2) — PDF merging and outline extraction

<br />
<br />
