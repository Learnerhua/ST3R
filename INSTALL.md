# Installation Guide

> Detailed step-by-step instructions for installing all dependencies of the ST3R pipeline.
> If you only need a quick overview, see the [README](README.md).
>
> **中文版安装指南**: [INSTALL.zh-CN.md](INSTALL.zh-CN.md)

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Install mamba (recommended)](#2-install-mamba-recommended)
3. [Create the 3 conda environments](#3-create-the-3-conda-environments)
4. [Install Chromium for Step 11](#4-install-chromium-for-step-11)
5. [GPU setup (optional)](#5-gpu-setup-optional)
6. [Verify installation](#6-verify-installation)
7. [Run scripts via conda wrapper](#7-run-scripts-via-conda-wrapper)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Prerequisites

| Component | Requirement | Notes |
|-----------|-------------|-------|
| OS | Linux (Ubuntu 20.04+ recommended) | macOS and WSL2 also work; Windows native is untested |
| Conda | Anaconda / Miniconda / Miniforge ≥23.x | [Miniforge](https://github.com/conda-forge/miniforge) is the recommended base distribution |
| Disk | ≥ 15 GB free | 3 envs combined ≈ 9-12 GB |
| RAM | ≥ 16 GB | For datasets with > 100K cells |
| GPU (optional) | NVIDIA GPU with CUDA 11.8+ driver | Steps 06 and 07 are GPU-accelerated |

---

## 2. Install mamba (recommended)

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

---

## 3. Create the 3 conda environments

The pipeline uses three conda environments, each tuned to the Python version its dependencies require:

| Env name | Python | Steps | Key packages |
|----------|--------|-------|--------------|
| `stereopy` | 3.8 | 01 | stereo, pandas |
| `scanpy` | 3.12 | 02-05, 11 | scanpy, squidpy, anndata, harmonypy, scrublet, jinja2, playwright, pypdf2 |
| `spateo_env` | 3.10 | 06-10 | spateo, torch (CUDA), pyvista |

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

---

## 4. Install Chromium for Step 11

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

---

## 5. GPU setup (optional)

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

## 6. Verify installation

Run these quick checks:

```bash
# Step 01 environment
conda run -n stereopy   python -c "import stereo, pandas; print('stereopy   : OK')"

# Steps 02-05 + 11 environment
conda run -n scanpy     python -c "import scanpy, squidpy, anndata, harmonypy, scrublet, jinja2; print('scanpy     : OK')"

# Steps 06-10 environment
conda run -n spateo_env python -c "import spateo, torch, pyvista, scanpy, anndata; print('spateo_env : OK | CUDA =', torch.cuda.is_available())"
```

All three should print `OK`.

---

## 7. Run scripts via conda wrapper

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
./run.sh 03_preprocess.py -I Output/02_concat/Sol_concated.h5ad -BK slice_id -P Sol_ -O Output/03_preprocess -minG 200 -maxG 2000 -minU 200 -maxU 6000 -minC 3 -maxMT 5 -maxHB 5
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

## 8. Troubleshooting

### `CondaHTTPError` / SSL errors
Configure a mirror (see step 2) or check your network proxy:
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

---

## Next steps

Return to the [README](README.md#usage) for step-by-step pipeline commands and the [Spateo_pipeline_SOP.pdf](Spateo_pipeline_SOP.pdf) for full parameter documentation.