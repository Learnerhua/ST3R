#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/python
import sys, argparse, warnings, platform, os, re, time, shutil
from datetime import datetime
import yaml


# ------------------------------------------------------------------ #
# Collect result files (PNG / interactive HTML / CSV / JSON)         #
# ------------------------------------------------------------------ #
def collect_results(root, prefix, sample_list_path, max_embed_mb, embed_enabled):
    """Locate upstream outputs under <root>/<NN>_* and return sections/slices/warnings."""
    sys.path.insert(0, _lib_dir())
    import collect
    return collect.collect(root, prefix, sample_list_path,
                           max_embed_mb, embed_enabled)


# ------------------------------------------------------------------ #
# Detect software versions used by upstream Scripts (01..10)         #
# ------------------------------------------------------------------ #
def detect_software_versions(scripts_dir, manual_versions):
    """Return a list of {name, version, desc} for the 6 main packages."""
    sys.path.insert(0, _lib_dir())
    import collect, content
    return collect.read_software_versions(
        scripts_dir, manual=manual_versions, desc_map=content.PACKAGE_DESCRIPTIONS)


# ------------------------------------------------------------------ #
# Render HTML (interactive + PDF-input) under -O                       #
# ------------------------------------------------------------------ #
def render_html_report(sections, slices, project, report_cfg, out_dir,
                       template_dir, assets_dir, software_versions,
                       workflow_svg, spateo_logo):
    """Build assets/, images/, viz/ under <out_dir> and render report.html + report_pdf.html."""
    sys.path.insert(0, _lib_dir())
    import render_html
    return render_html.render(
        sections, slices, project, report_cfg,
        out_dir, template_dir, assets_dir,
        software_versions=software_versions,
        workflow_svg=workflow_svg,
        spateo_logo=spateo_logo)


# ------------------------------------------------------------------ #
# Copy csv / json data files into <out_dir>/data/                     #
# ------------------------------------------------------------------ #
def copy_data_files(out_dir, sections):
    """Copy glm_data.csv and morph.json into <out_dir>/data/."""
    data_dir = os.path.join(out_dir, "data")
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    n_data = 0
    for sec in sections:
        for src in (sec.get("de_src"), sec.get("morph_src")):
            if src and os.path.isfile(src):
                shutil.copy2(src, os.path.join(data_dir, os.path.basename(src)))
                n_data += 1

    return data_dir, n_data


# ------------------------------------------------------------------ #
# Convert report_pdf.html -> report.pdf (Playwright)                  #
# ------------------------------------------------------------------ #
def generate_pdf(pdf_html, pdf_out, config, sections, cover_pdf, header_logo):
    """Render PDF with cover/TOC/header/footer/bookmarks via Playwright + PyPDF2."""
    sys.path.insert(0, _lib_dir())
    import html_to_pdf
    html_to_pdf.convert(pdf_html, pdf_out, config, sections,
                        cover_pdf, header_logo)


# ------------------------------------------------------------------ #
# Clean intermediate files (e.g. report_pdf.html)                     #
# ------------------------------------------------------------------ #
def cleanup_intermediates(out_dir, names):
    """Delete Playwright / build intermediates under <out_dir>."""
    removed = []
    for name in names:
        p = os.path.join(out_dir, name)
        if os.path.isfile(p):
            os.remove(p)
            removed.append(name)
    return removed


# ------------------------------------------------------------------ #
# Helpers                                                             #
# ------------------------------------------------------------------ #
def _lib_dir():
    """Return the Report_config/lib directory (sibling of Scripts/)."""
    return os.path.normpath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Report_config", "lib"))


def _bundle_paths(config_path):
    """Derive template / assets / cover / header paths from the config file path."""
    bundle = os.path.dirname(config_path)
    template_dir = os.path.join(bundle, "templates")
    assets_dir = os.path.join(bundle, "assets")
    cover_pdf = os.path.join(template_dir, "pdf_cover.pdf")
    header_logo = os.path.join(assets_dir, "images", "logo.png")
    return bundle, template_dir, assets_dir, cover_pdf, header_logo


def _print_run_info(root, prefix, sample_list_path, config_path, out_dir, no_pdf):
    """Print the standard run-information banner."""
    print("\n=========================== Report Generation ============================\n")
    print("Current time   :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Operating system:", platform.system(), platform.release())
    print("Platform       :", platform.platform())
    print("Working directory:", os.getcwd())
    print("Input root     :", root)
    print("File prefix    :", f"'{prefix}'" if prefix else "(none)")
    print("Sample list    :", sample_list_path or "(none)")
    print("Config file    :", config_path)
    print("Output directory:", out_dir)
    print("PDF generation :", "disabled (--no-pdf)" if no_pdf else "enabled")
    print("\nCopyright (c) 2026 KMHD. All Rights Reserved.")
    print("\n===========================================================================\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate the final analysis report (HTML + PDF) from the outputs of "
                    "scripts 01..10 of the Spateo 3D reconstruction pipeline. "
                    "This is step 11 of the pipeline; no upstream analysis is re-run."
    )

    parser.add_argument(
        "-I", "--input", required=True,
        metavar="INPUT_DIR",
        help="Root directory that contains the upstream output folders 01_* ... 10_*."
    )
    parser.add_argument(
        "-P", "--prefix", default="",
        metavar="PREFIX",
        help="Filename prefix shared by upstream results (default: none)."
    )
    parser.add_argument(
        "-SL", "--sample_list", default="",
        metavar="TSV_FILE",
        help="Optional TSV file listing sample_id / slice_id for each slice."
    )
    parser.add_argument(
        "-C", "--config", default="",
        metavar="YAML_FILE",
        help="Path to the report config.yaml (default: ../Report_config/config.yaml)."
    )
    parser.add_argument(
        "-O", "--outpath", default="./11_report",
        metavar="OUTPUT_DIR",
        help="Output directory for the report (default: ./11_report)."
    )
    parser.add_argument(
        "--no-pdf", action="store_true",
        help="Skip PDF generation; produce HTML only."
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    root              = os.path.abspath(args.input)
    prefix            = args.prefix
    sample_list_path  = os.path.abspath(args.sample_list) if args.sample_list else ""
    script_dir        = os.path.dirname(os.path.abspath(__file__))
    config_path       = os.path.abspath(args.config) if args.config else os.path.normpath(
        os.path.join(script_dir, "..", "Report_config", "config.yaml"))
    out_dir           = os.path.abspath(args.outpath)
    no_pdf            = args.no_pdf

    if not os.path.isdir(root):
        print(f"Error: Input directory not found: {root}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(config_path):
        print(f"Error: Config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)

    _, template_dir, assets_dir, cover_pdf, header_logo = _bundle_paths(config_path)

    # ------------------------------------------------------------------ #
    # Load config                                                          #
    # ------------------------------------------------------------------ #
    with open(config_path, "r", encoding="utf-8") as fh:
        config = yaml.safe_load(fh)
    project         = config.get("project", {}) or {}
    report_cfg      = config.get("report", {}) or {}
    embed_enabled   = bool(report_cfg.get("embed_interactive_3d", True))
    max_embed_mb    = float(report_cfg.get("max_embed_mb", 60))
    manual_versions = config.get("software_versions") or {}

    _print_run_info(root, prefix, sample_list_path, config_path, out_dir, no_pdf)

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Step 1: Locate upstream result files                                 #
    # ------------------------------------------------------------------ #
    print("[STEP 1] Locating upstream result files ...")
    sections, slices, warnings = collect_results(
        root, prefix, sample_list_path, max_embed_mb, embed_enabled)
    if not sections:
        print("Error: No displayable results found. Check -I / -P.", file=sys.stderr)
        sys.exit(1)
    n_fig = sum(len(s["figures"]) for s in sections)
    print(f"  -> {len(sections)} sections, {n_fig} figures, {len(slices)} slices")
    for w in warnings:
        print(f"  {w}")

    # ------------------------------------------------------------------ #
    # Step 1.5: Detect software versions                                   #
    # ------------------------------------------------------------------ #
    print("\n[STEP 1.5] Detecting software versions ...")
    software_versions = detect_software_versions(script_dir, manual_versions)
    n_manual = sum(1 for k, v in manual_versions.items() if v)
    print(f"  -> {len(software_versions)} entries "
          f"({n_manual} manual, {len(software_versions) - n_manual} auto)")

    # ------------------------------------------------------------------ #
    # Step 2: Render HTML (assets/, images/, viz/ under -O)                #
    # ------------------------------------------------------------------ #
    print("\n[STEP 2] Rendering HTML report ...")
    if os.path.exists(out_dir):
        for name in ("assets", "images", "viz", "data",
                     "report.html", "report.pdf"):
            p = os.path.join(out_dir, name)
            if os.path.isdir(p):
                shutil.rmtree(p)
            elif os.path.isfile(p):
                os.remove(p)
    else:
        os.makedirs(out_dir)

    workflow_svg = os.path.join(assets_dir, "images", "spateo_tdr_workflow.svg")
    spateo_logo  = os.path.join(assets_dir, "images", "spateo_logo.png")
    rendered = render_html_report(
        sections, slices, project, report_cfg, out_dir,
        template_dir, assets_dir, software_versions,
        workflow_svg, spateo_logo)
    print(f"  -> {rendered['html']}")

    # ------------------------------------------------------------------ #
    # Step 2.5: Copy csv / json data files                                 #
    # ------------------------------------------------------------------ #
    print("\n[STEP 2.5] Copying data files (csv / json) ...")
    data_dir, n_data = copy_data_files(out_dir, sections)
    print(f"  -> {data_dir} ({n_data} files)")

    # ------------------------------------------------------------------ #
    # Step 3: Generate PDF                                                 #
    # ------------------------------------------------------------------ #
    if not no_pdf:
        print("\n[STEP 3] Generating PDF report ...")
        try:
            pdf_out = os.path.join(out_dir, "report.pdf")
            generate_pdf(rendered["pdf_html"], pdf_out, config, sections,
                         cover_pdf, header_logo)
            size_mb = os.path.getsize(pdf_out) / (1024 * 1024)
            print(f"  -> {pdf_out} ({size_mb:.1f} MB)")
        except Exception as e:
            print(f"  [WARN] PDF generation failed (HTML report is intact): {e}")
            import traceback
            traceback.print_exc()
    else:
        print("\n[STEP 3] Skipped (--no-pdf)")

    # ------------------------------------------------------------------ #
    # Step 4: Per-section image subfolders                                 #
    #         (already produced in STEP 2: images/3.x_<anchor>/)           #
    # ------------------------------------------------------------------ #

    # ------------------------------------------------------------------ #
    # Step 5: Clean intermediate files                                     #
    # ------------------------------------------------------------------ #
    print("\n[STEP 5] Cleaning intermediate files ...")
    removed = cleanup_intermediates(out_dir, ("report_pdf.html",))
    for name in removed:
        print(f"  - removed: {name}")

    elapsed = time.time() - pipeline_start
    print("\n===========================================================================")
    print("Report generation complete!")
    print("Output directory:", out_dir)
    print("    report.html / report.pdf / images/ / viz/ / data/ / assets/")
    print("Total time   : {:.1f} s".format(elapsed))
    print("===========================================================================\n")


if __name__ == "__main__":
    main()