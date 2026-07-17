# -*- coding: utf-8 -*-
"""
content_subset.py
=================
Subset variant of ``content.py`` for the ``11_report_subset.py`` workflow.

Use case: input is a pre-annotated H5AD with cell-type labels (rather than the
full 01-10 pipeline output with squidpy-derived spatial domains). Sections 06
and 07 are rewritten to use "cell type" wording instead of "spatial domain",
because that is what the figures actually display in this workflow.

Mechanism:
    The lib modules (``collect.py``, ``render_html.py``) do ``import content``
    internally and reference ``content.RESULT_SECTIONS``. This module imports
    ``content`` and overwrites ``content.RESULT_SECTIONS`` with the
    cell-type wording, so those lib calls transparently pick up the override
    within the same Python process. Because ``content_subset`` is only
    imported by ``11_report_subset.py``, the full ``11_report.py`` (which
    never imports this module) keeps the original ``空间域`` wording.

All other names (METHODS / FAQ / REFERENCES / etc.) are inherited unchanged.
"""

# Re-export every other symbol from content so consumers can
# ``from content_subset import METHODS, FAQ, ...`` interchangeably.
from content import (  # noqa: F401  -- re-export
    SPATEO_INTRO,
    METHODS_INTRO,
    METHODS_WORKFLOW_SVG,
    SPATEO_LOGO,
    METHODS,
    PACKAGE_DESCRIPTIONS,
    MORPH_LABELS,
    DELIVERY,
    FAQ,
    REFERENCES,
)

# ------------------------------------------------------------------ #
# Subset-specific RESULT_SECTIONS — only 06 and 07 are rewritten.     #
# ------------------------------------------------------------------ #
RESULT_SECTIONS = [
    # ----------------------------- 06 对齐 -----------------------------
    {
        "key": "06",
        "anchor": "alignment",
        "title": "切片空间对齐",
        "intro": (
            "时空转录组连续切片在切片、染色与成像过程中会产生平移、旋转甚至"
            "局部形变。若直接按 Z 轴堆叠，三维模型会出现明显的“错位”与“扭曲”，"
            "后续基于此的三维重建将不可靠。本项目输入为已经完成细胞类型注释的"
            "H5AD 数据，每个 spot 自带细胞类型标签；Spateo 通过生成式高斯过程"
            "与最优传输算法，在保持基因表达一致性与细胞类型空间分布一致性的"
            "前提下，将相邻切片配准到统一坐标系，是构建高精度三维模型的关键"
            "前提。"
        ),
        "figures": [
            {"suffix": "squidpy_2Dslices", "title": "对齐前各切片空间分布",
             "caption": "每个子图为一张切片，按原始空间坐标绘制，颜色表示细胞类型；用于展示对齐前切片间的位置关系。"},
            {"suffix": "aligned_2Dslices", "title": "对齐后各切片空间分布",
             "caption": "坐标含义同上，展示经空间对齐后各切片被配准到统一坐标系的结果。"},
            {"suffix": "aligned_2Dslices_overlap", "title": "相邻切片对齐效果（随机 3 张示例）",
             "caption": "随机抽取 3 张相邻切片叠加于同一坐标系，不同颜色代表不同切片，用于直观展示对齐后相邻切片的空间重合程度与形变校正效果。"},
        ],
    },
    # ----------------------------- 07 三维重建 -----------------------------
    {
        "key": "07",
        "anchor": "tdr",
        "title": "三维组织模型重建",
        "intro": (
            "在二维分析只能刻画每一切片“平面上的表达结构”的前提下，三维重建进一步"
            "沿 Z 轴堆叠对齐后的切片，恢复组织在深度方向的真实形态。本项目基于 Spateo "
            "依次重建三类模型：(1) 三维点云（point cloud）模型——每个 spot 还原为三维"
            "空间中的一个点，颜色可映射至其细胞类型或表达量，可直观呈现各细胞类型"
            "在组织三维空间中的分布；(2) 三维网格（mesh）曲面模型——在点云基础上重建"
            "组织外表面，便于直观观察组织整体形态；(3) 三维体素（voxel）模型——将三维"
            "空间离散化为规则体素网格，可用于体积统计与内部结构分区展示。"
        ),
        "figures": [
            {"suffix": "aligned_pc_3D", "title": "三维点云模型", "interactive": True,
             "caption": "三维散点为 spot，坐标为重建后的三维空间坐标，颜色表示细胞类型；可交互旋转观察组织整体三维结构。"},
            {"suffix": "aligned_pc_3D_multi", "title": "三维点云模型（多视角）",
             "caption": "从多个视角展示同一三维点云模型，颜色表示细胞类型，便于从不同方向观察组织形态。"},
            {"suffix": "aligned_pc_3D_multi_clusters", "title": "三维点云模型（各细胞类型分视图）",
             "caption": "每个子图单独展示一种细胞类型在三维空间中的分布，用于观察各细胞类型的三维位置与形态。"},
            {"suffix": "aligned_mesh_3D", "title": "三维网格曲面模型", "interactive": True,
             "caption": "由点云重建得到的组织外表面网格（mesh）模型，展示组织的三维轮廓与体表形态。"},
            {"suffix": "aligned_voxel_3D", "title": "三维体素模型", "interactive": True,
             "caption": "将三维模型体素化（voxel）后的实体表示，颜色表示细胞类型，用于展示组织内部的三维分区结构。"},
        ],
    },
]


# ------------------------------------------------------------------ #
# Override content.RESULT_SECTIONS for lib modules (collect / render) #
# ------------------------------------------------------------------ #
import content as _content
_content.RESULT_SECTIONS = RESULT_SECTIONS