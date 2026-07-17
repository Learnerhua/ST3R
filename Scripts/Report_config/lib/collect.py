# -*- coding: utf-8 -*-
"""
collect.py
==========
确定性地定位前 10 步流程在结果根目录下产出的文件，并组织成供模板渲染的数据结构。

不依赖智能体：全部按 "结果根目录 / {两位数编号}_* / {前缀}{后缀}.png" 的命名规则匹配；
文件缺失时自动跳过并记录警告，保证任意项目、任意前缀、任意切片数均可运行。
"""

import os
import re
import glob
import json
import csv
import subprocess

import content


# -----------------------------------------------------------------------------
# 软件版本探测
# -----------------------------------------------------------------------------
# 11 步运行时不重跑分析，但仍需要在报告里展示本次实际使用的软件版本。
# 思路：读 Scripts/01..10.py 的 shebang，按其指向的 conda 环境启动该环境下的 python，
# 探测指定包的 __version__。任何探测失败都不抛错，调用方按需降级为"未知"。
#
# 报告只展示以下 6 个主要软件：
_DISPLAY_PACKAGES = [
    ("stereo",   "Stereopy"),
    ("scanpy",   "Scanpy"),
    ("squidpy",  "Squidpy"),
    ("spateo",   "Spateo"),
    ("torch",    "PyTorch"),
    ("pyvista",  "PyVista"),
]

# 各脚本里会用到的包，用于按 shebang 自动发现其环境
_SCRIPT_PKG_HINTS = {
    "01_gef2h5ad.py":     ["stereo"],
    "02_concat.py":       ["anndata"],
    "03_preprocess.py":   ["scanpy"],
    "04_squidpy.py":      ["squidpy", "scanpy"],
    "05_dataConvert.py":  ["scanpy"],
    "06_align.py":        ["spateo", "torch", "scanpy"],
    "07_tdr.py":          ["spateo", "pyvista", "torch"],
    "08_backbone.py":     ["spateo", "pyvista"],
    "09_morph.py":        ["spateo", "pyvista"],
    "10_interpolation.py": ["spateo", "scanpy", "pyvista"],
}


def _parse_shebang_env(shebang_line):
    """从 shebang 行解析 conda 环境名。
    形如 '#!/oldhome/ouyjh/miniforge3/envs/Stereopy/bin/python' -> 'Stereopy'。
    """
    m = re.search(r"/envs/([^/]+)/bin/python", shebang_line or "")
    return m.group(1) if m else None


def _detect_version(env_python, package_name):
    """在该环境 python 中 import 包并取 __version__。失败返回 None。"""
    if not env_python or not os.path.exists(env_python):
        return None
    try:
        out = subprocess.check_output(
            [env_python, "-c",
             f"import {package_name} as p; print(getattr(p, '__version__', 'unknown'))"],
            stderr=subprocess.DEVNULL,
            timeout=15,
        )
        return out.decode().strip() or None
    except (subprocess.SubprocessError, OSError):
        return None


def _resolve_env_for_pkg(scripts_dir, conda_root, pkg_import):
    """在 Scripts/01..10.py 中找使用 pkg_import 的脚本，返回对应环境名。"""
    for fname, pkgs in _SCRIPT_PKG_HINTS.items():
        if pkg_import not in pkgs:
            continue
        path = os.path.join(scripts_dir, fname)
        if not os.path.exists(path):
            continue
        try:
            with open(path, "r", encoding="utf-8") as fh:
                first = fh.readline()
        except OSError:
            continue
        env = _parse_shebang_env(first)
        if env:
            return env
    return None


def read_software_versions(scripts_dir, conda_root="/oldhome/ouyjh/miniforge3/envs",
                           manual=None, desc_map=None):
    """
    收集 6 个主要软件的版本。优先级：
        manual[display_name]   客户在 config.yaml 里手动填的（最优先）
        否则自动按 shebang 探测（探测失败返回 "未知"）
        都没有则不展示该软件

    Returns
    -------
    list of dict: [{'name': 'Spateo', 'version': '0.x.y', 'desc': '...'}, ...]
    """
    manual = manual or {}
    desc_map = desc_map or {}
    out = []
    for pkg_import, display_name in _DISPLAY_PACKAGES:
        version = None
        if display_name in manual and manual[display_name]:
            version = str(manual[display_name]).strip()
        else:
            env = _resolve_env_for_pkg(scripts_dir, conda_root, pkg_import)
            if env:
                env_python = os.path.join(conda_root, env, "bin", "python")
                version = _detect_version(env_python, pkg_import)
        if version:
            out.append({"name": display_name, "version": version,
                        "desc": desc_map.get(display_name, "")})
    return out


def _find_step_dir(root, key):
    """在结果根目录下按 '{key}_*' 匹配某一步骤的输出子目录。"""
    matches = sorted(glob.glob(os.path.join(root, f"{key}_*")))
    matches = [m for m in matches if os.path.isdir(m)]
    return matches[0] if matches else None


def _file_size_mb(path):
    try:
        return os.path.getsize(path) / (1024 * 1024)
    except OSError:
        return 0.0


def parse_sample_list(sample_list_path):
    """
    解析 sample_list.tsv，返回按 slice_id 自然排序的 [(slice_id, sample_id), ...]。
    需要列包含 sample_id 与 slice_id（gef_path 可选）。
    """
    if not sample_list_path or not os.path.exists(sample_list_path):
        return []

    rows = []
    with open(sample_list_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            sid = (r.get("sample_id") or "").strip()
            slid = (r.get("slice_id") or "").strip()
            if sid or slid:
                rows.append((slid, sid))

    def _natural_key(item):
        slid = item[0]
        # 提取末尾数字用于自然排序（Sol_1, Sol_2, ... Sol_10）
        num = ""
        for ch in reversed(slid):
            if ch.isdigit():
                num = ch + num
            else:
                break
        return (slid[: len(slid) - len(num)], int(num) if num else 0, slid)

    rows.sort(key=_natural_key)
    return rows


def _resolve_figure(step_dir, prefix, fig, max_embed_mb, embed_enabled):
    """解析单张图：定位 png 与可选的交互式 html，判断是否内嵌。"""
    suffix = fig["suffix"]
    png = os.path.join(step_dir, f"{prefix}{suffix}.png")
    if not os.path.exists(png):
        return None

    item = {
        "title": fig.get("title", suffix),
        "caption": fig.get("caption", ""),
        "png_path": png,
        "png_name": os.path.basename(png),
        "interactive_path": None,   # 内嵌用的交互式 html 绝对路径
    }

    if fig.get("interactive"):
        html = os.path.join(step_dir, f"{prefix}{suffix}.html")
        if os.path.exists(html):
            size = _file_size_mb(html)
            if embed_enabled and size <= max_embed_mb:
                item["interactive_path"] = html
    return item


def _load_morph(step_dir, prefix, morph_suffix):
    """读取形态学指标 json，返回 [(中文名, 值, 释义), ...]。"""
    path = os.path.join(step_dir, f"{prefix}{morph_suffix}.json")
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return None

    rows = []
    for k, v in data.items():
        label, desc = content.MORPH_LABELS.get(k, (k, ""))
        if isinstance(v, float):
            av = abs(v)
            if av == 0:
                vstr = "0"
            elif av < 1e-3:
                # 极小浮点数（如 cell_density = 4e-05）保留有效数字，避免被四舍五入到 0
                vstr = f"{v:.3g}"
            elif av < 1e6:
                vstr = f"{v:,.4f}".rstrip("0").rstrip(".") or "0"
            else:
                vstr = f"{v:,.2f}"
        else:
            vstr = str(v)
        rows.append((label, vstr, desc))
    return rows


def _load_de_table(step_dir, csv_name, top_n=20):
    """读取差异表达 csv，返回 (表头, 前 top_n 行, 总行数)。"""
    path = os.path.join(step_dir, csv_name)
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = next(reader, [])
            rows = [row for row in reader]
    except OSError:
        return None
    return {"header": header, "rows": rows[:top_n], "total": len(rows), "name": csv_name}


def collect(root, prefix, sample_list_path, max_embed_mb, embed_enabled):
    """
    汇总全部结果，返回 (sections, slices, warnings)。
    sections : 供模板渲染的章节列表（仅包含实际存在图片的章节）。
    slices   : [(slice_id, sample_id), ...]
    warnings : 缺失文件的告警字符串列表。
    """
    warnings = []
    sections = []

    for sec in content.RESULT_SECTIONS:
        step_dir = _find_step_dir(root, sec["key"])
        if step_dir is None:
            warnings.append(f"[跳过] 未找到步骤 {sec['key']} 的输出目录（{root}/{sec['key']}_*）")
            continue

        figures = []
        for fig in sec.get("figures", []):
            resolved = _resolve_figure(step_dir, prefix, fig, max_embed_mb, embed_enabled)
            if resolved is None:
                warnings.append(
                    f"[缺图] {sec['title']}：未找到 {prefix}{fig['suffix']}.png"
                )
            else:
                figures.append(resolved)

        morph = None
        morph_src = None
        if sec.get("morph_json"):
            morph_src = os.path.join(step_dir, f"{prefix}{sec['morph_json']}.json")
            if os.path.exists(morph_src):
                morph = _load_morph(step_dir, prefix, sec["morph_json"])
            else:
                warnings.append(f"[缺表] {sec['title']}：未找到形态学指标 json")
                morph_src = None

        de_table = None
        de_src = None
        if sec.get("de_csv"):
            de_src = os.path.join(step_dir, sec["de_csv"])
            if os.path.exists(de_src):
                de_table = _load_de_table(step_dir, sec["de_csv"])
            else:
                warnings.append(f"[缺表] {sec['title']}：未找到差异表达结果 {sec['de_csv']}")
                de_src = None

        # 该章节没有任何可展示内容则跳过
        if not figures and not morph and not de_table:
            warnings.append(f"[跳过] {sec['title']}：无可展示的结果文件")
            continue

        sections.append({
            "key": sec["key"],
            "anchor": sec["anchor"],
            "title": sec["title"],
            "intro": sec["intro"],
            "note": sec.get("note", ""),
            "figures": figures,
            "morph": morph,
            "morph_src": morph_src,
            "de_table": de_table,
            "de_src": de_src,
        })

    slices = parse_sample_list(sample_list_path)
    if not slices:
        warnings.append("[提示] 未提供或未解析到切片信息（sample_list.tsv）")

    return sections, slices, warnings
