# 安装指南

> ST3R 流水线所有依赖的详细安装步骤。
> 如只需快速概览，请参见 [README.zh-CN.md](README.zh-CN.md)（[English README](README.md)）。

---

## 目录

1. [前置条件](#1-前置条件)
2. [安装 mamba（推荐）](#2-安装-mamba推荐)
3. [创建 3 个 conda 环境](#3-创建-3-个-conda-环境)
4. [为 Step 11 安装 Chromium](#4-为-step-11-安装-chromium)
5. [GPU 配置（可选）](#5-gpu-配置可选)
6. [验证安装](#6-验证安装)
7. [通过 conda 包装脚本运行脚本](#7-通过-conda-包装脚本运行脚本)
8. [故障排查](#8-故障排查)

---

## 1. 前置条件

| 组件 | 要求 | 备注 |
|------|------|------|
| 操作系统 | Linux（Ubuntu 20.04+ 推荐） | macOS 与 WSL2 也支持；Windows 原生未测试 |
| Conda | Anaconda / Miniconda / Miniforge ≥ 23.x | 推荐使用 [Miniforge](https://github.com/conda-forge/miniforge) 作为基础发行版 |
| 磁盘 | ≥ 15 GB 可用空间 | 3 个环境合计约 9-12 GB |
| 内存 | ≥ 16 GB | 用于 > 100K 细胞的数据集 |
| GPU（可选） | 支持 CUDA 11.8+ 的 NVIDIA GPU | Steps 06 和 07 支持 GPU 加速 |

---

## 2. 安装 mamba（推荐）

`mamba` 是 `conda` 的 C++ 重实现，依赖解析速度快 10 倍。在基础环境中安装：

```bash
conda install -n base -c conda-forge mamba -y
```

若已使用 [Miniforge](https://github.com/conda-forge/miniforge) 或 [Mambaforge](https://github.com/conda-forge/miniforge)，`mamba` 已默认可用。

国内用户配置清华镜像源加速下载：

```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
conda config --set show_channel_urls yes
```

---

## 3. 创建 3 个 conda 环境

流水线使用三个 conda 环境，每个环境根据其依赖所要求的 Python 版本进行调优：

| 环境名 | Python | 步骤 | 关键包 |
|--------|--------|------|--------|
| `stereopy` | 3.8 | 01 | stereo, pandas |
| `scanpy` | 3.12 | 02-05, 11 | scanpy, squidpy, anndata, harmonypy, scrublet, jinja2, playwright, pypdf2 |
| `spateo_env` | 3.10 | 06-10 | spateo, torch (CUDA), pyvista |

在仓库根目录下依次执行：

```bash
# 1. Stereopy 环境（Step 01：GEF → H5AD）
mamba env create -f envs/stereopy.yml

# 2. Scanpy 环境（Steps 02-05 + Step 11：预处理、空间域、报告）
mamba env create -f envs/scanpy.yml

# 3. Spateo 环境（Steps 06-10：3D 重建核心）
mamba env create -f envs/spateo_env.yml
```

每个命令耗时 5-15 分钟，取决于网速。

若倾向使用 `conda`（较慢但内置）：

```bash
conda env create -f envs/stereopy.yml
conda env create -f envs/scanpy.yml
conda env create -f envs/spateo_env.yml
```

---

## 4. 为 Step 11 安装 Chromium

Step 11 通过 Playwright 调用 headless Chromium 将 HTML 渲染为 PDF。创建 `scanpy` 环境后，下载浏览器二进制：

```bash
conda run -n scanpy playwright install chromium
```

服务器（Ubuntu/Debian）可选安装 headless Chromium 所需的系统库：

```bash
sudo apt install -y libnss3 libnspr4 libatk1.0-0 libatk-bridge2.0-0 \
    libcups2 libdrm2 libxkbcommon0 libxcomposite1 libxdamage1 \
    libxfixes3 libxrandr2 libgbm1 libpango-1.0-0 libcairo2 libasound2
```

---

## 5. GPU 配置（可选）

Steps 06（`06_align.py`）和 07（`07_tdr.py`）使用 CUDA 加速。`spateo_env.yml` 中已包含来自 `pytorch` channel 的 `pytorch::pytorch`。

若您拥有 NVIDIA GPU：

1. 安装 NVIDIA 驱动（CUDA 12.x 需 ≥ 525；CUDA 11.8 需 ≥ 520）。
2. 验证 CUDA 对 PyTorch 可见：
   ```bash
   conda run -n spateo_env python -c "import torch; print('CUDA:', torch.cuda.is_available(), '|', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'no GPU')"
   ```
   有 GPU 时的预期输出：`CUDA: True | NVIDIA GeForce RTX ...`
3. 若**没有 GPU**，编辑 `envs/spateo_env.yml`，移除 `pytorch` channel 与 `pytorch::pytorch` 行，改用 PyPI 上的 CPU 版 torch：
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

## 6. 验证安装

运行以下快速检查：

```bash
# Step 01 环境
conda run -n stereopy   python -c "import stereo, pandas; print('stereopy   : OK')"

# Steps 02-05 + 11 环境
conda run -n scanpy     python -c "import scanpy, squidpy, anndata, harmonypy, scrublet, jinja2; print('scanpy     : OK')"

# Steps 06-10 环境
conda run -n spateo_env python -c "import spateo, torch, pyvista, scanpy, anndata; print('spateo_env : OK | CUDA =', torch.cuda.is_available())"
```

三项均应输出 `OK`。

---

## 7. 通过 conda 包装脚本运行脚本

随仓库分发的脚本使用**硬编码的 shebang**，指向占位符路径（如 `#!/path/to/envs/Stereopy/bin/python`）。您必须将每个 shebang 改为本地路径，或使用随仓库根目录一同分发的 `run.sh` 包装脚本。

仓库根目录已自带 `run.sh` 包装脚本，通过正确的 conda 环境调用每个脚本。您无需自行创建，只需赋予执行权限：

```bash
chmod +x run.sh
```

（脚本源码如下供参考，您无需复制。）

```bash
#!/bin/bash
# run.sh — 使用正确的 conda 环境调用 ST3R 脚本。

set -e

if [ $# -lt 1 ]; then
    echo "Usage: ./run.sh <script_name.py> [args...]"
    echo "Example: ./run.sh 01_gef2h5ad.py -C rawData/sample.tsv -BT cell_bins -O Output/01"
    exit 1
fi

SCRIPT="$1"
shift

# 映射脚本名 → conda 环境
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
        echo "[run.sh] 未知脚本: $SCRIPT"
        exit 1
        ;;
esac

# 通过 conda run 调用脚本，并将仓库的 Scripts/ 目录加入 PYTHONPATH
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec conda run -n "$ENV" \
    env PYTHONPATH="$SCRIPT_DIR/Scripts/Report_config/lib:${PYTHONPATH:-}" \
    python "$SCRIPT_DIR/Scripts/$SCRIPT" "$@"
```

使用示例：

```bash
./run.sh 01_gef2h5ad.py -C rawData/sample_list.tsv -BT cell_bins -O Output/01_gef2h5ad
./run.sh 03_preprocess.py -I Output/02_concat/Sol_concated.h5ad -BK slice_id -P Sol_ -O Output/03_preprocess -minG 200 -maxG 2000 -minU 200 -maxU 6000 -minC 3 -maxMT 5 -maxHB 5
./run.sh 07_tdr.py -AD Output/06_alignment/Sol_adata_aligned.h5ad -RD Output/05_dataConvert/Sol_compatible.h5ad -P Sol_ -O Output/07_tdr
./run.sh 11_report.py -I Output -P Sol_ -SL rawData/sample_list.tsv -O Output/11_report
```

> **为什么不直接修改 shebang？** 编辑 12 个脚本中的 shebang 会让每个安装耦合到硬编码的 Python 路径。`run.sh` 包装脚本保持环境布局的灵活性。

### 备选方案：直接修改 shebang 行

如果不想使用包装脚本——例如希望直接通过 `bash Scripts/01_gef2h5ad.py` 调用脚本——可以直接编辑每个脚本的 shebang 行，使其指向本地的 Python 路径。每个脚本只需编辑一次，对应其所属的 conda 环境：

| 脚本 | 替换 shebang 为 |
|------|----------------|
| `01_gef2h5ad.py` | `#!/path/to/your/envs/stereopy/bin/python` |
| `02_concat.py`、`03_preprocess.py`、`04_squidpy.py`、`05_dataConvert.py` | `#!/path/to/your/envs/scanpy/bin/python` |
| `06_align.py`、`07_tdr.py`、`08_backbone.py`、`09_morph.py`、`10_interpolation.py` | `#!/path/to/your/envs/spateo_env/bin/python` |
| `11_report.py`、`11_report_subset.py` | `#!/path/to/your/envs/scanpy/bin/python` |

通过 `conda env list` 或 `mamba env list` 查找路径，输出形如：

```
stereopy                 /path/to/envs/stereopy
scanpy                   /path/to/envs/scanpy
spateo_env               /path/to/envs/spateo_env
```

然后通过一条 sed 命令批量替换所有脚本：

```bash
# 将旧 shebang 替换为本地路径
sed -i '1c\#!/path/to/envs/stereopy/bin/python'   Scripts/01_gef2h5ad.py
sed -i '1c\#!/path/to/envs/scanpy/bin/python'     Scripts/02_concat.py Scripts/03_preprocess.py Scripts/04_squidpy.py Scripts/05_dataConvert.py Scripts/11_report.py Scripts/11_report_subset.py
sed -i '1c\#!/path/to/envs/spateo_env/bin/python' Scripts/06_align.py Scripts/07_tdr.py Scripts/08_backbone.py Scripts/09_morph.py Scripts/10_interpolation.py
```

之后即可直接调用脚本：

```bash
./Scripts/01_gef2h5ad.py -C rawData/sample_list.tsv -BT cell_bins -O Output/01_gef2h5ad
./Scripts/07_tdr.py -AD Output/06_alignment/Sol_adata_aligned.h5ad -RD Output/05_dataConvert/Sol_compatible.h5ad -P Sol_ -O Output/07_tdr
```

> **权衡**：修改 shebang 在运行时更简洁（无需包装脚本），但脚本会与特定机器耦合。重新克隆或共享仓库时会覆盖您的修改。若希望仓库保持可移植性，建议使用 `run.sh`。

---

## 8. 故障排查

### `CondaHTTPError` / SSL 错误
配置镜像源（见第 2 节）或检查网络代理：
```bash
conda config --show channels
mamba clean --all
```

### `PackageNotFoundError` / 依赖求解失败
最常见的原因是缺少 channel。确保 `.condarc` 包含：
```yaml
channels:
  - conda-forge
  - bioconda
  - pytorch
show_channel_urls: true
```

### `pytorch::pytorch` 安装失败（spateo_env）
尝试从 PyPI 安装 torch：
```bash
conda activate spateo_env
pip install --extra-index-url https://download.pytorch.org/whl/cu118 torch torchvision
pip install spateo
```

### `playwright install chromium` 失败 / Chromium 无法启动
1. 使用 `--with-deps` 重新运行（Linux 上需 `sudo`）：
   ```bash
   conda run -n scanpy playwright install --with-deps chromium
   ```
2. 若处于公司防火墙后，运行前设置 `HTTPS_PROXY`。

### GPU 未被识别
1. 验证 NVIDIA 驱动：
   ```bash
   nvidia-smi
   ```
2. 确认 CUDA 工具包版本与 PyTorch 期望的版本匹配（如 `cu118` wheels 对应 CUDA 11.8）。
3. 从匹配的 `--extra-index-url` 重新安装 torch。

### 运行 `11_report.py` 时 ImportError
该脚本从 `Scripts/Report_config/lib/` 中导入 `content` 和 `content_subset`。上述 `run.sh` 包装脚本已正确设置 `PYTHONPATH`。若直接调用 Python，需手动指定：
```bash
export PYTHONPATH="$PWD/Scripts/Report_config/lib:$PYTHONPATH"
conda run -n scanpy python Scripts/11_report.py ...
```

### 需要重建环境
```bash
mamba env remove -n stereopy
mamba env create -f envs/stereopy.yml
```

---

## 下一步

返回 [README.zh-CN.md](README.zh-CN.md)（[English README](README.md)）查阅逐步的流水线命令，参见 [Spateo_pipeline_SOP.pdf](Spateo_pipeline_SOP.pdf) 获取完整参数文档。