# Installation Guide

> Detailed step-by-step instructions for installing all dependencies of the ST3R pipeline.
> If you only need a quick overview, see the [README](README.md).
>
> **中文版安装指南**: [INSTALL.zh-CN.md](INSTALL.zh-CN.md)

---

## Table of Contents

1. [Choose an installation method](#1-choose-an-installation-method)
2. [Method A — Use the prebuilt Docker / Singularity image (recommended)](#2-method-a--use-the-prebuilt-docker--singularity-image-recommended)
3. [Method B — Install conda environments manually](#3-method-b--install-conda-environments-manually)
4. [Verify installation](#4-verify-installation)
5. [Run scripts](#5-run-scripts)
6. [Troubleshooting](#6-troubleshooting)

---

## 1. Choose an installation method

The pipeline ships with **two equivalent installation paths**:

| Method | Best for | Time to first run |
|--------|----------|-------------------|
| **A. Prebuilt image** (Docker / Singularity) | Most users; cross-platform; reproducible | ~5 minutes (only image pull) |
| **B. Manual conda envs** from `envs/*.yml` | Users who can't run containers, or need to modify dependencies | ~30-45 minutes (full env build) |

If you have Docker or Singularity installed, **Method A is strongly recommended**. The image is identical to the environment used to produce the reference results shipped in `Sol_test/`, so all parameters have been validated against it.

---

## 2. Method A — Use the prebuilt Docker / Singularity image (recommended)

A prebuilt image with all 3 conda environments (`stereopy`, `scanpy`, `spateo_env`) is available on Docker Hub, and the corresponding Singularity `.sif` is archived on Zenodo for long-term reference:

> **Docker Hub** — `oyjhlovedocker/spateo_tdr:v1`
> <https://hub.docker.com/r/oyjhlovedocker/spateo_tdr>
>
> **Zenodo** — `spateo_tdr.sif` (DOI: [10.5281/zenodo.21776415](https://doi.org/10.5281/zenodo.21776415))
> <https://doi.org/10.5281/zenodo.21776415>

The image is built from this repo's `envs/*.yml` plus the `spateo-release` source tree, so it reproduces the exact environment used to validate the pipeline.

### 2.1 Pull the image

```bash
# Docker — pull directly from Docker Hub
docker pull oyjhlovedocker/spateo_tdr:v1

# Singularity / Apptainer — either convert on the fly from Docker Hub
singularity pull spateo_tdr.sif docker://oyjhlovedocker/spateo_tdr:v1

# ...or download the prebuilt .sif directly from Zenodo (recommended for HPC
# clusters without Docker Hub access)
wget -O spateo_tdr.sif https://zenodo.org/records/21776415/files/spateo_tdr.sif
```

### 2.2 Verify the image

```bash
# Docker: list envs inside the container
docker run --rm oyjhlovedocker/spateo_tdr:v1 bash -c \
    "conda env list && which python"

# Singularity
singularity exec --cleanenv --no-mount tmp spateo_tdr.sif conda env list
```

> **Why `--cleanenv` for every step?**
> - `--cleanenv` clears inherited host env vars (`PYTHONPATH`, `LD_LIBRARY_PATH`, etc.) so the container's conda envs activate cleanly.
>
> **Why `--no-mount tmp` only for GPU steps (06-10)?**
> - The `spateo_env` conda environment was installed into the container's `/tmp` directory at build time. If Singularity bind-mounts the host's `/tmp` over the container's `/tmp`, the conda env is masked and `conda run -n spateo_env ...` fails. `--no-mount tmp` is required **only when invoking `spateo_env`**, i.e. Steps 06-10.

Expected: three envs (`stereopy`, `scanpy`, `spateo_env`), each with its own Python.

### 2.3 Run a single pipeline step

The image ships only the pre-installed conda environments; it does **not** contain the pipeline repo, raw data, or output. You mount host paths at runtime with `-v` (Docker) or `-B` (Singularity). The examples below use these conventions (replace with your actual paths):

| Container path | Host path | Purpose |
|---------------|-----------|---------|
| `/Pipeline` | `/path/to/Pipeline` | Repo root (contains `Scripts/`, `Report_config/`, `envs/`, `run.sh`) |
| `/rawData` | `/path/to/rawData` | Input data (GEF files, `sample_list.tsv`) — mounted read-only |
| `/Output` | `/path/to/Output` | Output dir for all step artifacts — mounted read-write |

> **Prerequisite**: This section assumes familiarity with container volume mounts (`-v` / `-B`) and GPU passthrough (`--gpus all` / `--nv`). If you are new to Docker / Singularity, consult their official docs first.

```bash
# === Docker ===
docker run --rm \
    -v /path/to/Pipeline:/Pipeline \
    -v /path/to/rawData:/rawData:ro \
    -v /path/to/Output:/Output \
    oyjhlovedocker/spateo_tdr:v1 \
    conda run -n stereopy python /Pipeline/Scripts/01_gef2h5ad.py \
        -C /rawData/sample_list.tsv -BT cell_bins -O /Output/01_gef2h5ad

# === Singularity / Apptainer ===
singularity exec --cleanenv --no-mount tmp \
    -B /path/to/Pipeline:/Pipeline \
    -B /path/to/rawData:/rawData:ro \
    -B /path/to/Output:/Output \
    spateo_tdr.sif \
    conda run -n stereopy python /Pipeline/Scripts/01_gef2h5ad.py \
        -C /rawData/sample_list.tsv -BT cell_bins -O /Output/01_gef2h5ad
```

### 2.4 Run GPU-accelerated steps (06 / 07 / 08 / 09 / 10)

Steps 06-10 use CUDA. Pass the GPU through to the container:

```bash
# Docker — use NVIDIA Container Toolkit
docker run --rm --gpus all \
    -v /path/to/Pipeline:/Pipeline \
    -v /path/to/Output:/Output \
    oyjhlovedocker/spateo_tdr:v1 \
    conda run -n spateo_env python /Pipeline/Scripts/07_tdr.py \
        -AD /Output/06_alignment/Sol_adata_aligned.h5ad \
        -RD /Output/05_dataConvert/Sol_compatible.h5ad \
        -P Sol_ -O /Output/07_tdr

# Singularity — use --nv + --no-mount tmp (spateo_env installed in /tmp)
singularity exec --cleanenv --no-mount tmp --nv \
    -B /path/to/Pipeline:/Pipeline \
    -B /path/to/Output:/Output \
    spateo_tdr.sif \
    conda run -n spateo_env python /Pipeline/Scripts/07_tdr.py \
        -AD /Output/06_alignment/Sol_adata_aligned.h5ad \
        -RD /Output/05_dataConvert/Sol_compatible.h5ad \
        -P Sol_ -O /Output/07_tdr
```

### 2.5 Run the full pipeline in one shell

Run Steps 01-11 sequentially via a simple shell loop, activating the correct env per step. Example:

```bash
PIPELINE=/path/to/Pipeline
OUTPUT=/path/to/Output
SLICE_LIST=/path/to/rawData/sample_list.tsv

docker run --rm -v $PIPELINE:$PIPELINE -v /path/to/rawData:/rawData:ro -v $OUTPUT:$OUTPUT \
    oyjhlovedocker/spateo_tdr:v1 bash -c "
        conda run -n stereopy   python $PIPELINE/Scripts/01_gef2h5ad.py -C /rawData/sample_list.tsv -BT cell_bins -O $OUTPUT/01_gef2h5ad
        conda run -n scanpy     python $PIPELINE/Scripts/02_concat.py -I $OUTPUT/01_gef2h5ad -O $OUTPUT/02_concat -P Sol_ -minC 3 -maxMT 10
        conda run -n scanpy     python $PIPELINE/Scripts/03_preprocess.py -I $OUTPUT/02_concat/Sol_concated.h5ad -BK batch -P Sol_ -O $OUTPUT/03_preprocess
        conda run -n scanpy     python $PIPELINE/Scripts/04_squidpy.py -I $OUTPUT/03_preprocess/Sol_preprocessed.h5ad -O $OUTPUT/04_squidpy -LK slice_id -P Sol_ -R 1.2 -WS 0.8
        conda run -n scanpy     python $PIPELINE/Scripts/05_dataConvert.py -I $OUTPUT/04_squidpy/Sol_squidpy.h5ad -P Sol_ -O $OUTPUT/05_dataConvert
        conda run -n spateo_env python $PIPELINE/Scripts/06_align.py -I $OUTPUT/05_dataConvert/Sol_compatible.h5ad -P Sol_ -O $OUTPUT/06_alignment
        conda run -n spateo_env python $PIPELINE/Scripts/07_tdr.py -AD $OUTPUT/06_alignment/Sol_adata_aligned.h5ad -RD $OUTPUT/05_dataConvert/Sol_compatible.h5ad -P Sol_ -O $OUTPUT/07_tdr
        conda run -n spateo_env python $PIPELINE/Scripts/08_backbone.py -AD $OUTPUT/07_tdr/Sol_tdr.h5ad -PC $OUTPUT/07_tdr/Sol_aligned_pc_model.vtk -MS $OUTPUT/07_tdr/Sol_aligned_mesh_model.vtk -P Sol_ -O $OUTPUT/08_backbone
        conda run -n spateo_env python $PIPELINE/Scripts/09_morph.py -PC $OUTPUT/07_tdr/Sol_aligned_pc_model.vtk -MS $OUTPUT/07_tdr/Sol_aligned_mesh_model.vtk -O $OUTPUT/09_Morph -P Sol_
        conda run -n spateo_env python $PIPELINE/Scripts/10_interpolation.py -AD $OUTPUT/07_tdr/Sol_tdr.h5ad -PC $OUTPUT/07_tdr/Sol_aligned_pc_model.vtk -MS $OUTPUT/07_tdr/Sol_aligned_mesh_model.vtk -VX $OUTPUT/07_tdr/Sol_aligned_voxel_model.vtk -GL $OUTPUT/08_backbone/glm_data.csv -O $OUTPUT/10_interpolation -P Sol_ -NG 3 -NS 15
        conda run -n scanpy     python $PIPELINE/Scripts/11_report.py -I $OUTPUT -P Sol_ -SL /rawData/sample_list.tsv -O $OUTPUT/11_report
    "
```

(Singularity users: replace the `docker run ...` prefix with `singularity exec --cleanenv -B ... oyjh_spateo_tdr_v1.sif bash -c "..."`. For lines that run `conda run -n spateo_env ...` (Steps 06-10), also add `--no-mount tmp --nv` — see [§2.2](#22-verify-the-image) for why.)

---

## 3. Method B — Install conda environments manually

Use this method if you can't run containers or need to customize dependency versions.

### 3.1 Prerequisites

| Component | Requirement | Notes |
|-----------|-------------|-------|
| OS | Linux (Ubuntu 20.04+ recommended) | macOS and WSL2 also work; Windows native is untested |
| Conda | Anaconda / Miniconda / Miniforge ≥ 23.x | [Miniforge](https://github.com/conda-forge/miniforge) is the recommended base distribution |
| Disk | ≥ 15 GB free | 3 envs combined ≈ 9-12 GB |
| RAM | ≥ 16 GB | For datasets with > 100K cells |
| GPU (optional) | NVIDIA GPU with CUDA 11.8+ driver | Steps 06 and 07 are GPU-accelerated |

### 3.2 Install mamba (recommended)

`mamba` is a C++ reimplementation of `conda` and is 10× faster at solving dependencies. Install it into your base environment:

```bash
conda install -n base -c conda-forge mamba -y
```

If you already use [Miniforge](https://github.com/conda-forge/miniforge) or [Mambaforge](https://github.com/conda-forge/miniforge), `mamba` is already available.

Configure the Tsinghua mirror (China users) for faster downloads:

```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
conda config --set show_channel_urls yes
```

### 3.3 Create the 3 conda environments

The pipeline uses three conda environments, each tuned to the Python version its dependencies require:

| Env name | Python | Steps | Key packages |
|----------|--------|-------|--------------|
| `stereopy` | 3.8 | 01 | stereopy, pandas |
| `scanpy` | 3.12 | 02-05, 11 | scanpy, squidpy, anndata, jinja2, playwright, pypdf2 |
| `spateo_env` | 3.10 | 06-10 | spateo-release, torch (CUDA), pyvista |

Run these commands from the repository root:

```bash
# 1. Stereopy environment (Step 01: GEF → H5AD)
mamba env create -f envs/stereopy.yml

# 2. Scanpy environment (Steps 02-05 + Step 11: preprocessing, domains, reporting)
mamba env create -f envs/scanpy.yml

# 3. Spateo environment (Steps 06-10: 3D reconstruction core)
mamba env create -f envs/spateo_env.yml
```

Each command takes 5-15 minutes depending on network speed.

If you prefer `conda` (slower but built-in):

```bash
conda env create -f envs/stereopy.yml
conda env create -f envs/scanpy.yml
conda env create -f envs/spateo_env.yml
```

### 3.4 Install Chromium for Step 11

Step 11 renders HTML to PDF using headless Chromium via Playwright. After creating the `scanpy` environment, fetch the browser binary:

```bash
conda run -n scanpy playwright install chromium
```

Optional: install system libraries required by headless Chromium on a server (Ubuntu/Debian):

```bash
sudo apt install -y libnss3 libnspr4 libatk1.0-0 libatk-bridge2.0-0 \
    libcups2 libdrm2 libxkbcommon0 libxcomposite1 libxdamage1 \
    libxfixes3 libxrandr2 libgbm1 libpango-1.0-0 libcairo2 libasound2
```

### 3.5 GPU setup (optional)

Steps 06 (`06_align.py`) and 07 (`07_tdr.py`) use CUDA for acceleration. The `spateo_env.yml` already includes `pytorch::pytorch` from the `pytorch` channel.

If you have an NVIDIA GPU:

1. Install the NVIDIA driver (≥ 525 for CUDA 12.x; ≥ 520 for CUDA 11.8).
2. Verify CUDA is visible to PyTorch:
   ```bash
   conda run -n spateo_env python -c "import torch; print('CUDA:', torch.cuda.is_available(), '|', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'no GPU')"
   ```
   Expected output with a GPU: `CUDA: True | NVIDIA GeForce RTX ...`
3. If you have **no GPU**, edit `envs/spateo_env.yml` and remove the `pytorch` channel + the `pytorch::pytorch` lines, then rely on the CPU-only torch from PyPI:
   ```yaml
   dependencies:
     - python=3.10
     - pip
     - pip:
         - --extra-index-url https://download.pytorch.org/whl/cpu
         - torch
         - torchvision
         - spateo
   ```

---

## 4. Verify installation

Run these quick checks:

**For Method A (container image):**

```bash
docker run --rm oyjhlovedocker/spateo_tdr:v1 bash -c '
    conda run -n stereopy   python -c "import stereopy, pandas; print(\"stereopy   : OK\")"
    conda run -n scanpy     python -c "import scanpy, squidpy, anndata, jinja2; print(\"scanpy     : OK\")"
    conda run -n spateo_env python -c "import spateo, torch, pyvista; print(\"spateo_env : OK | CUDA =\", torch.cuda.is_available())"
'
```

**For Method B (manual conda install):**

```bash
# Step 01 environment
conda run -n stereopy   python -c "import stereopy, pandas; print('stereopy   : OK')"

# Steps 02-05 + 11 environment
conda run -n scanpy     python -c "import scanpy, squidpy, anndata, jinja2; print('scanpy     : OK')"

# Steps 06-10 environment
conda run -n spateo_env python -c "import spateo, torch, pyvista, scanpy, anndata; print('spateo_env : OK | CUDA =', torch.cuda.is_available())"
```

All three should print `OK`.

---

## 5. Run scripts

### 5.1 If you used Method A (container)

Invoke scripts directly via `docker run` or `singularity exec` (see [§2.3](#23-run-a-single-pipeline-step) and [§2.4](#24-run-gpu-accelerated-steps-06--07--08--09--10)). The image ships with all 3 conda envs pre-activated by name (`stereopy`, `scanpy`, `spateo_env`).

> Singularity users: always pass `--cleanenv` to every `singularity exec`. Pass `--no-mount tmp --nv` **only** to Steps 06-10 (the `spateo_env` was installed in `/tmp` at build time — see [§2.2](#22-verify-the-image)).

For a full one-shot driver covering all 11 steps, see [§2.5](#25-run-the-full-pipeline-in-one-shell).

### 5.2 If you used Method B (manual conda)

The shipped scripts use **hardcoded shebangs** that point to a placeholder path (e.g. `#!/path/to/envs/Stereopy/bin/python`). You must either edit each shebang to your local path, or use the bundled `run.sh` wrapper (already shipped in the repository root).

A `run.sh` wrapper ships in the repository root that invokes each script via the correct conda environment. Just make it executable:

```bash
chmod +x run.sh
```

The shipped scripts use **hardcoded shebangs** that point to a placeholder path (e.g. `#!/path/to/envs/Stereopy/bin/python`). You must either edit each shebang to your local path, or use the bundled `run.sh` wrapper (already shipped in the repository root).

A `run.sh` wrapper ships in the repository root that invokes each script via the correct conda environment. Just make it executable:

```bash
chmod +x run.sh
```

For reference, here is the wrapper's source (already shipped — no need to copy):

```bash
#!/bin/bash
# run.sh — Invoke ST3R scripts using the correct conda environment.

set -e

if [ $# -lt 1 ]; then
    echo "Usage: ./run.sh <script_name.py> [args...]"
    echo "Example: ./run.sh 01_gef2h5ad.py -C rawData/sample.tsv -BT cell_bins -O Output/01"
    exit 1
fi

SCRIPT="$1"
shift

# Map script name → conda environment
case "$(basename "$SCRIPT" .py)" in
    01_gef2h5ad)
        ENV=stereopy
        ;;
    02_concat|03_preprocess|04_squidpy|05_dataConvert|11_report|11_report_subset)
        ENV=scanpy
        ;;
    06_align|07_tdr|08_backbone|09_morph|10_interpolation)
        ENV=spateo_env
        ;;
    *)
        echo "[run.sh] Unknown script: $SCRIPT"
        exit 1
        ;;
esac

# Invoke script via conda run, with the repo's Scripts/ dir on PYTHONPATH
# (needed for Step 11 to find Report_config/lib/*).
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec conda run -n "$ENV" \
    env PYTHONPATH="$SCRIPT_DIR/Scripts/Report_config/lib:${PYTHONPATH:-}" \
    python "$SCRIPT_DIR/Scripts/$SCRIPT" "$@"
```

```bash
./run.sh 01_gef2h5ad.py -C rawData/sample_list.tsv -BT cell_bins -O Output/01_gef2h5ad
./run.sh 03_preprocess.py -I Output/02_concat/Sol_concated.h5ad -BK batch -P Sol_ -O Output/03_preprocess
./run.sh 07_tdr.py -AD Output/06_alignment/Sol_adata_aligned.h5ad -RD Output/05_dataConvert/Sol_compatible.h5ad -P Sol_ -O Output/07_tdr
./run.sh 11_report.py -I Output -P Sol_ -SL rawData/sample_list.tsv -O Output/11_report
```

> **Why not just fix the shebangs?** Editing the shebangs in 12 script files couples every install to a hard-coded Python path. The `run.sh` wrapper keeps your environment layout flexible.

### Alternative: Fix the shebang lines directly

If you prefer not to use a wrapper script — for example, to invoke scripts directly via `bash Scripts/01_gef2h5ad.py` — you can edit each script's shebang to point to your local Python path. Each script must be edited once to use the matching environment's interpreter:

| Script | Replace shebang with |
|--------|----------------------|
| `01_gef2h5ad.py` | `#!/path/to/your/envs/stereopy/bin/python` |
| `02_concat.py`, `03_preprocess.py`, `04_squidpy.py`, `05_dataConvert.py` | `#!/path/to/your/envs/scanpy/bin/python` |
| `06_align.py`, `07_tdr.py`, `08_backbone.py`, `09_morph.py`, `10_interpolation.py` | `#!/path/to/your/envs/spateo_env/bin/python` |
| `11_report.py`, `11_report_subset.py` | `#!/path/to/your/envs/scanpy/bin/python` |

Find the path with `conda env list` or `mamba env list`, which prints something like:

```
stereopy                 /path/to/envs/stereopy
scanpy                   /path/to/envs/scanpy
spateo_env               /path/to/envs/spateo_env
```

Then run a one-shot sed across all scripts:

```bash
# Replace the old shebangs with your local paths
sed -i '1c\#!/path/to/envs/stereopy/bin/python'   Scripts/01_gef2h5ad.py
sed -i '1c\#!/path/to/envs/scanpy/bin/python'     Scripts/02_concat.py Scripts/03_preprocess.py Scripts/04_squidpy.py Scripts/05_dataConvert.py Scripts/11_report.py Scripts/11_report_subset.py
sed -i '1c\#!/path/to/envs/spateo_env/bin/python' Scripts/06_align.py Scripts/07_tdr.py Scripts/08_backbone.py Scripts/09_morph.py Scripts/10_interpolation.py
```

After this, you can invoke scripts directly:

```bash
./Scripts/01_gef2h5ad.py -C rawData/sample_list.tsv -BT cell_bins -O Output/01_gef2h5ad
./Scripts/07_tdr.py -AD Output/06_alignment/Sol_adata_aligned.h5ad -RD Output/05_dataConvert/Sol_compatible.h5ad -P Sol_ -O Output/07_tdr
```

> **Trade-off**: Editing shebangs is simpler at runtime (no wrapper) but couples the scripts to a specific machine. Re-cloning or sharing the repo will overwrite your edits. Use `run.sh` if you want the repo to remain portable.

---

## 6. Troubleshooting

### `CondaHTTPError` / SSL errors
Configure a mirror (see [§3.2](#32-install-mamba-recommended)) or check your network proxy:
```bash
conda config --show channels
mamba clean --all
```

### `PackageNotFoundError` / dependency solver fails
The most common cause is missing channels. Ensure your `.condarc` includes:
```yaml
channels:
  - conda-forge
  - bioconda
  - pytorch
show_channel_urls: true
```

### `pytorch::pytorch` install fails (spateo_env)
Try installing torch from PyPI instead:
```bash
conda activate spateo_env
pip install --extra-index-url https://download.pytorch.org/whl/cu118 torch torchvision
pip install spateo
```

### `playwright install chromium` fails / Chromium won't launch
1. Re-run with `--with-deps` (requires `sudo` on Linux):
   ```bash
   conda run -n scanpy playwright install --with-deps chromium
   ```
2. If behind a corporate firewall, set `HTTPS_PROXY` before running.

### GPU not detected
1. Verify the NVIDIA driver:
   ```bash
   nvidia-smi
   ```
2. Ensure CUDA toolkit matches PyTorch's expected version (e.g., `11.8` for `cu118` wheels).
3. Reinstall torch from the matching `--extra-index-url`.

### ImportError when running `11_report.py`
The script imports `content` and `content_subset` from `Scripts/Report_config/lib/`. The `run.sh` wrapper above sets `PYTHONPATH` correctly. If you invoke Python directly, prepend it manually:
```bash
export PYTHONPATH="$PWD/Scripts/Report_config/lib:$PYTHONPATH"
conda run -n scanpy python Scripts/11_report.py ...
```

### Need to recreate an environment
```bash
mamba env remove -n stereopy
mamba env create -f envs/stereopy.yml
```

### Docker: `permission denied` while binding mounts
On Linux, add your user to the `docker` group, or prefix the command with `sudo`. Files written into mounted volumes will be owned by `root` unless you also pass `--user "$(id -u):$(id -g)"`.

### Singularity: `FATAL: container creation failed`
If the host kernel is too old or missing required capabilities, try `--writable-tmpfs` or update Singularity to ≥ 3.10.

### GPU not visible inside the container
- **Docker**: install [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) and restart the Docker daemon. Verify with `docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi`.
- **Singularity**: ensure the host has the NVIDIA driver loaded (`nvidia-smi` works on host), then use `singularity exec --cleanenv --no-mount tmp --nv ...` for GPU steps (06-10). Steps 01-05 and 11 only need `--cleanenv`.

---

## Next steps

Return to the [README](README.md#usage) for step-by-step pipeline commands and the [Spateo_pipeline_SOP.pdf](Spateo_pipeline_SOP.pdf) for full parameter documentation.