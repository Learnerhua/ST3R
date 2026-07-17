# -*- coding: utf-8 -*-
"""
render_html.py
==============
将定位到的结果数据渲染为可浏览的多文件 HTML（报告构建目录）。

产出（在 build 目录下）：
    assets/        report.css, logo.png, KMHD_logo_circle.png
    images/        全部结果图片（png）
    viz/           需内嵌的交互式三维模型（html）
    report.html    交互版（三维模型以 iframe 展示）
    report_pdf.html 打印版（三维模型以静态截图展示，供转 PDF 用）
"""

import os
import shutil
from jinja2 import Environment, FileSystemLoader

import content


def _copy_unique(src, dst_dir):
    """复制文件到目标目录，返回文件名（若重名则保留原名，假定同名即同文件）。"""
    name = os.path.basename(src)
    dst = os.path.join(dst_dir, name)
    if not os.path.exists(dst):
        shutil.copy2(src, dst)
    return name


def render(sections, slices, project, report_cfg, build_dir, template_dir, assets_dir,
           software_versions=None, spateo_intro=None, workflow_svg=None, spateo_logo=None):
    build_assets = os.path.join(build_dir, "assets")
    build_images = os.path.join(build_dir, "images")
    build_viz = os.path.join(build_dir, "viz")
    for d in (build_assets, build_images, build_viz):
        os.makedirs(d, exist_ok=True)

    # 复制静态资源
    shutil.copy2(os.path.join(assets_dir, "css", "report.css"),
                 os.path.join(build_assets, "report.css"))
    for img in ("logo.png", "KMHD_logo_circle.png"):
        src = os.path.join(assets_dir, "images", img)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(build_assets, img))

    # 复制图片与交互式模型，并回填文件名（含按步骤子目录的相对路径）
    for sec_idx, sec in enumerate(sections, start=1):
        sub = os.path.join(build_images, f"3.{sec_idx}_{sec['anchor']}")
        os.makedirs(sub, exist_ok=True)
        sub_rel = f"images/3.{sec_idx}_{sec['anchor']}"
        for fig in sec["figures"]:
            png_basename = _copy_unique(fig["png_path"], sub)
            fig["png_name"] = png_basename            # 文件名（仅 basename）
            fig["png_relpath"] = f"{sub_rel}/{png_basename}"  # HTML 引用的相对路径
            if fig.get("interactive_path"):
                fig["viz_name"] = _copy_unique(fig["interactive_path"], build_viz)

    # 复制 Spateo logo / workflow SVG（缺资源时优雅跳过）
    workflow_svg_exists = False
    if workflow_svg and os.path.exists(workflow_svg):
        shutil.copy2(workflow_svg, os.path.join(build_assets, os.path.basename(workflow_svg)))
        workflow_svg_rel = "assets/" + os.path.basename(workflow_svg)
        workflow_svg_exists = True
    else:
        workflow_svg_rel = ""

    spateo_logo_exists = False
    if spateo_logo and os.path.exists(spateo_logo):
        shutil.copy2(spateo_logo, os.path.join(build_assets, os.path.basename(spateo_logo)))
        spateo_logo_rel = "assets/" + os.path.basename(spateo_logo)
        spateo_logo_exists = True
    else:
        spateo_logo_rel = ""

    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=False,   # 图注等为可信静态内容
    )
    template = env.get_template("report.html.j2")

    common = dict(
        project=project,
        report=report_cfg,
        slices=slices,
        sections=sections,
        methods_intro=content.METHODS_INTRO,
        methods=content.METHODS,
        spateo_intro=spateo_intro or content.SPATEO_INTRO,
        workflow_svg=workflow_svg_rel,
        workflow_svg_exists=workflow_svg_exists,
        spateo_logo=spateo_logo_rel,
        spateo_logo_exists=spateo_logo_exists,
        software_versions=software_versions or [],
        delivery=content.DELIVERY,
        faq=content.FAQ,
        references=content.REFERENCES,
    )

    html_path = os.path.join(build_dir, "report.html")
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(template.render(mode="html", **common))

    pdf_html_path = os.path.join(build_dir, "report_pdf.html")
    with open(pdf_html_path, "w", encoding="utf-8") as fh:
        fh.write(template.render(mode="pdf", **common))

    return {"html": html_path, "pdf_html": pdf_html_path,
            "images_dir": build_images, "viz_dir": build_viz}
