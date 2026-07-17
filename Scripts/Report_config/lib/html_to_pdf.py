# -*- coding: utf-8 -*-
"""
html_to_pdf.py
==============
使用 Playwright 将打印版 HTML（report_pdf.html，三维模型为静态图）转换为 PDF，
并叠加封页、目录、页眉页脚与 PDF 书签大纲。

改编自 KMHD FASQR 报告流程，目录页由实际章节结构动态生成。
"""

import os
import re
import base64
from typing import List, Dict


def _image_to_base64(image_path):
    if not image_path or not os.path.exists(image_path):
        return ""
    with open(image_path, "rb") as f:
        enc = base64.b64encode(f.read()).decode("utf-8")
    ext = os.path.splitext(image_path)[1].lower()
    mime = {"png": "png", "jpg": "jpeg", "jpeg": "jpeg", "svg": "svg+xml"}.get(
        ext.lstrip("."), "png")
    return f"data:image/{mime};base64,{enc}"


# ------------------------- PDF 大纲提取 -------------------------
def _extract_outline(pdf_path: str):
    try:
        from PyPDF2 import PdfReader
        reader = PdfReader(pdf_path)

        def process(items, level=0):
            out = []
            i = 0
            while i < len(items):
                item = items[i]
                if isinstance(item, dict) and "/Title" in item:
                    try:
                        page = reader.get_destination_page_number(item) + 1
                    except Exception:
                        page = 1
                    node = {"title": str(item["/Title"]).strip(), "page": page,
                            "level": level, "children": []}
                    if i + 1 < len(items) and isinstance(items[i + 1], list):
                        node["children"] = process(items[i + 1], level + 1)
                        i += 1
                    out.append(node)
                elif isinstance(item, list):
                    out.extend(process(item, level))
                i += 1
            return out

        if hasattr(reader, "outline") and reader.outline:
            return process(reader.outline)
    except Exception as e:
        print(f"[WARN] 提取 PDF 大纲失败: {e}")
    return []


def _title_page_map(tree) -> Dict[str, int]:
    m = {}

    def walk(items):
        for it in items:
            m.setdefault(it["title"], it["page"])
            if it.get("children"):
                walk(it["children"])
    walk(tree)
    return m


# ------------------------- 目录页 -------------------------
def _build_toc_entries(sections):
    """根据实际章节结构生成 [(标题, 层级), ...]。"""
    entries = [("1. 项目信息", 1), ("2. 分析方法", 1), ("3. 分析结果", 1)]
    for i, sec in enumerate(sections, start=1):
        entries.append((f"3.{i} {sec['title']}", 2))
    entries += [
        ("4. 报告解读与交付说明", 1),
        ("4.1 交付文件说明", 2),
        ("4.2 交付范围说明", 2),
        ("5. 常见问题", 1),
        ("6. 参考文献", 1),
    ]
    return entries


def _create_toc_html(entries, title_to_page):
    items = []
    for title, level in entries:
        page = title_to_page.get(title, "")
        colors = {1: ("#4682B4", "#2c3e50", "600"),
                  2: ("#5D9CEC", "#34495e", "500"),
                  3: ("#7FB3D5", "#5D6D7E", "normal")}
        ncol, tcol, fw = colors.get(level, colors[3])
        mt = re.match(r"^(\d+(?:\.\d+)*)[.\s]*(.+)", title)
        if mt:
            colored = (f'<span class="toc-number" style="color:{ncol};font-weight:{fw};">'
                       f'{mt.group(1)}.</span>'
                       f'<span class="toc-text" style="color:{tcol};">{mt.group(2)}</span>')
        else:
            colored = f'<span class="toc-text" style="color:{tcol};font-weight:{fw};">{title}</span>'
        indent = (level - 1) * 20
        items.append(f'''
        <li class="toc-item level-{level}">
          <div class="toc-line">
            <div class="toc-title-container" style="margin-left:{indent}px;">{colored}</div>
            <span class="toc-dots"></span>
            <span class="toc-page" style="color:#4682B4;font-weight:500;">{page}</span>
          </div>
        </li>''')
    items_html = "\n".join(items)
    return f'''<!DOCTYPE html><html><head><meta charset="UTF-8"><style>
    @page {{ margin: 15mm 20mm; }}
    body {{ font-family:"Noto Sans CJK SC","Microsoft YaHei",Arial,sans-serif; margin:0; color:#333; }}
    .toc-header {{ text-align:center; margin-bottom:10mm; padding-bottom:5mm; border-bottom:2px solid #4682B4; }}
    .toc-title {{ font-size:22pt; font-weight:bold; color:#4682B4; margin:0; letter-spacing:3px; }}
    .toc-subtitle {{ font-size:11pt; color:#7FB3D5; margin-top:3mm; font-style:italic; }}
    .toc-list {{ list-style:none; padding:0; margin:5mm 0 0 0; }}
    .toc-item {{ margin:8px 0; line-height:1.5; }}
    .toc-line {{ display:flex; justify-content:space-between; align-items:baseline; }}
    .toc-title-container {{ display:flex; align-items:baseline; background:#fff; padding-right:8px; }}
    .toc-number {{ margin-right:6px; min-width:25px; }}
    .level-1 .toc-text {{ font-size:11.5pt; }}
    .level-2 .toc-text {{ font-size:10.5pt; }}
    .toc-dots {{ flex:1; border-bottom:1px dotted #B0C4DE; margin:0 10px; top:-3px; position:relative; opacity:.6; }}
    .toc-page {{ background:#fff; padding-left:8px; min-width:18px; text-align:right; }}
    </style></head><body>
    <div class="toc-header"><h1 class="toc-title">目 录</h1><div class="toc-subtitle">Table of Contents</div></div>
    <ul class="toc-list">{items_html}</ul></body></html>'''


# ------------------------- 页眉页脚 -------------------------
def _header_template(cfg, logo_b64):
    if not cfg.get("enabled", True):
        return "<div></div>"
    logo = f'<img src="{logo_b64}" style="height:{cfg.get("logo_height","25px")};" />' if logo_b64 else ""
    right = cfg.get("right_text", "")
    return f'''<div style="width:100%;font-size:10px;padding:5px 0 0 0;display:flex;
        align-items:flex-start;justify-content:space-between;box-sizing:border-box;position:relative;">
        <div style="margin-left:10mm;">{logo}</div>
        <div style="color:#333;font-size:9px;text-align:right;margin-right:10mm;align-self:flex-end;">{right}</div>
        <div style="position:absolute;bottom:0;left:10mm;right:10mm;border-top:1px solid #000;height:0;"></div>
        </div>'''


def _footer_template(cfg):
    if not cfg.get("enabled", True):
        return "<div></div>"
    text = cfg.get("text", "")
    pages = ""
    if cfg.get("show_page_numbers", True):
        pages = ('<div style="position:absolute;left:50%;transform:translateX(-50%);'
                 'color:#4682B4;font-weight:500;white-space:nowrap;">'
                 '第 <span class="pageNumber"></span> 页 / 共 <span class="totalPages"></span> 页</div>')
    return f'''<div style="width:100%;font-size:9px;padding:8px 0 5px 0;display:flex;
        align-items:center;justify-content:space-between;color:#666;box-sizing:border-box;position:relative;">
        <div style="margin-left:10mm;flex:1;">{text}</div>{pages}
        <div style="margin-right:10mm;flex:1;"></div></div>'''


# ------------------------- PDF 合并 -------------------------
def _merge_pdfs(cover, main, out, toc):
    from PyPDF2 import PdfReader, PdfWriter
    writer = PdfWriter()
    if cover and os.path.exists(cover):
        writer.append_pages_from_reader(PdfReader(cover))
    if toc and os.path.exists(toc):
        writer.append_pages_from_reader(PdfReader(toc))
    main_start = len(writer.pages)
    for page in PdfReader(main).pages:
        writer.add_page(page)

    tree = _extract_outline(main)

    def add_tree(items, parent=None):
        for it in items:
            node = writer.add_outline_item(
                title=it["title"], page_number=it["page"] - 1 + main_start, parent=parent)
            if it.get("children"):
                add_tree(it["children"], node)
    if tree:
        add_tree(tree)

    with open(out, "wb") as fh:
        writer.write(fh)


# ------------------------- 主入口 -------------------------
def convert(pdf_html, output_pdf, config, sections, cover_pdf, header_logo):
    from playwright.sync_api import sync_playwright

    pdf_cfg = config.get("pdf_header_footer", {})
    header_cfg = pdf_cfg.get("header", {})
    footer_cfg = pdf_cfg.get("footer", {})
    margin_cfg = pdf_cfg.get("margins", {})
    margins = {"top": margin_cfg.get("top", "25mm"), "bottom": margin_cfg.get("bottom", "20mm"),
               "left": margin_cfg.get("left", "10mm"), "right": margin_cfg.get("right", "10mm")}

    logo_b64 = _image_to_base64(header_logo) if header_cfg.get("enabled", True) else ""
    header_html = _header_template(header_cfg, logo_b64)
    footer_html = _footer_template(footer_cfg)

    out_dir = os.path.dirname(os.path.abspath(output_pdf))
    temp_main = os.path.join(out_dir, "_temp_main.pdf")
    temp_toc = os.path.join(out_dir, "_temp_toc.pdf")

    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True, args=["--disable-dev-shm-usage", "--disable-gpu", "--no-sandbox"])

        # 主报告
        page = browser.new_page()
        page.goto(f"file://{os.path.abspath(pdf_html)}", wait_until="networkidle")
        page.wait_for_timeout(800)
        page.pdf(path=temp_main, format="A4", margin=margins, print_background=True,
                 prefer_css_page_size=True, tagged=True, outline=True,
                 display_header_footer=True, header_template=header_html,
                 footer_template=footer_html)
        page.close()

        # 目录页（据主报告大纲页码生成）
        tree = _extract_outline(temp_main)
        title_to_page = _title_page_map(tree)
        entries = _build_toc_entries(sections)
        toc_html = _create_toc_html(entries, title_to_page)
        tpage = browser.new_page()
        tpage.set_content(toc_html, wait_until="networkidle")
        tpage.wait_for_timeout(500)
        tpage.pdf(path=temp_toc, format="A4", margin=margins, print_background=True,
                  display_header_footer=False)
        tpage.close()
        browser.close()

    cover = cover_pdf if (cover_pdf and os.path.exists(cover_pdf)) else None
    _merge_pdfs(cover, temp_main, output_pdf, temp_toc)

    for t in (temp_main, temp_toc):
        if os.path.exists(t):
            os.remove(t)
    return output_pdf
