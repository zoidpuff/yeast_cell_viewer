"""Render a Markdown document (with local images) to PDF through headless Chromium.

  python scripts/md_to_pdf.py IN.md OUT.pdf

Images are resolved relative to the Markdown file. Needs the `markdown` package and Playwright's Chromium; on the
cluster set LD_LIBRARY_PATH/FONTCONFIG_FILE for .chromelibs (HANDOVER, headless browser notes).
"""
import sys, os, asyncio, pathlib
import markdown

CSS = """
@page { size: A4; margin: 18mm 16mm 18mm 16mm; }
body { font-family: 'DejaVu Sans', 'Liberation Sans', sans-serif; font-size: 10pt; line-height: 1.45; color: #1a1a19; }
h1 { font-size: 20pt; margin: 0 0 6pt; }
h2 { font-size: 14pt; margin: 18pt 0 6pt; padding-bottom: 3pt; border-bottom: 1px solid #d8d7d2; page-break-after: avoid; }
h3 { font-size: 11.5pt; margin: 12pt 0 4pt; page-break-after: avoid; }
p, li { orphans: 3; widows: 3; }
code { font-family: 'DejaVu Sans Mono', monospace; font-size: 8.5pt; background: #f0efec; padding: 0 2px; border-radius: 2px; }
table { border-collapse: collapse; margin: 6pt 0 10pt; font-size: 8.8pt; }
tr { page-break-inside: avoid; }
th, td { border-bottom: 1px solid #e3e2dc; padding: 3pt 6pt; text-align: left; vertical-align: top; }
th { border-bottom: 1.5px solid #a3a29c; font-weight: 600; }
img { max-width: 100%; display: block; margin: 6pt auto 10pt; page-break-inside: avoid; }
li > table { margin-left: 0; }
"""


async def main(src, dst):
    src = pathlib.Path(src).resolve()
    html = markdown.markdown(src.read_text(), extensions=["tables", "fenced_code", "sane_lists"])
    doc = f"<!doctype html><html><head><meta charset='utf-8'><base href='{src.parent.as_uri()}/'><style>{CSS}</style></head><body>{html}</body></html>"
    tmp = src.parent / (".render_" + src.stem + ".html")
    tmp.write_text(doc)
    from playwright.async_api import async_playwright
    async with async_playwright() as p:
        b = await p.chromium.launch(headless=True)
        pg = await b.new_page()
        await pg.goto(tmp.as_uri()); await pg.wait_for_load_state("networkidle")
        await pg.pdf(path=dst, format="A4", print_background=True, display_header_footer=True,
                     header_template="<span></span>",
                     footer_template="<div style='font-size:7pt;color:#777;width:100%;text-align:center'><span class='pageNumber'></span> / <span class='totalPages'></span></div>",
                     margin=dict(top="16mm", bottom="18mm", left="16mm", right="16mm"))
        await b.close()
    tmp.unlink()


if __name__ == "__main__":
    asyncio.run(main(sys.argv[1], sys.argv[2]))
