"""Simulated GFP micrographs of marker proteins (IMPROVEMENTS I2): a midplane gallery through the nucleus and maximum
projections of whole-cell z-stacks, plus foci counts. Optics as in pipeline/microscopy.py (widefield, NA 1.4 oil).
Output: work/microscopy/gallery_midplane.png, mip_punctate.png, spot_counts.json
"""
import os, sys, json
import numpy as np
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK
import microscopy as m

OUT = f"{WORK}/microscopy"
GALLERY = [("TDH3", "cytosol"), ("NOP1", "nucleolus"), ("NUP49", "nuclear pore"), ("HTB2", "chromatin"),
           ("PMA1", "plasma membrane"), ("PIL1", "eisosome"), ("SEC61", "ER"), ("ATP2", "mitochondria"),
           ("VPH1", "vacuole membrane"), ("VMA2", "V-ATPase V1 (vacuole face)"), ("POT1", "peroxisome"),
           ("SEC21", "COPI coat (Golgi face)"), ("SEC7", "Golgi"), ("ERG6", "lipid droplet"), ("ABP1", "actin patch"),
           ("DCP2", "P-body (log phase)")]
MIP = ["POT1", "SEC7", "SEC21", "PIL1", "SUR7", "ERG6", "ABP1", "DCP2", "SPC42", "VMA2"]


def tile(img, label, sub):
    a = np.asarray(img, float); lo, hi = np.percentile(a, 1), np.percentile(a, 99.8)
    u8 = (np.clip((a - lo) / max(hi - lo, 1e-9), 0, 1) ** 0.8 * 255).astype(np.uint8).T[::-1]
    I = Image.fromarray(u8).resize((u8.shape[1] * 3, u8.shape[0] * 3), Image.NEAREST).convert("RGB")
    I = Image.merge("RGB", (Image.new("L", I.size, 0), I.split()[1], Image.new("L", I.size, 0)))
    d = ImageDraw.Draw(I); d.text((4, 2), label, fill=(255, 255, 255)); d.text((4, 14), sub, fill=(200, 200, 200))
    d.line([(I.size[0] - 50, I.size[1] - 8), (I.size[0] - 4, I.size[1] - 8)], fill=(255, 255, 255), width=3)  # 1 um
    return I


def sheet(tiles, cols, path):
    W, H = tiles[0].size; rows = (len(tiles) + cols - 1) // cols
    S = Image.new("RGB", (cols * W + (cols - 1) * 4, rows * H + (rows - 1) * 4), (40, 40, 40))
    for i, t in enumerate(tiles):
        S.paste(t, ((i % cols) * (W + 4), (i // cols) * (H + 4)))
    S.save(path)


def main():
    os.makedirs(OUT, exist_ok=True)
    c = m.Cell(); nuc = c.pos[c.region == 1].mean(0)
    tiles = []
    for g, lab in GALLERY:
        p = c.positions(g)
        tiles.append(tile(m.render(c, p, z=nuc[2], photons=40, seed=1), f"{g} ({len(p)})", lab))
    sheet(tiles, 4, f"{OUT}/gallery_midplane.png")
    st, _ = m.render_stack(c, c.positions("TDH3")); ref = st.max(0); mask = ref > 0.25 * np.percentile(ref, 99)
    tiles, counts = [], {}
    for g in MIP:
        s, _ = m.render_stack(c, c.positions(g), photons=40, seed=2)
        counts[g] = m.count_spots(s, mask)[0]
        tiles.append(tile(s.max(0), f"{g} MIP", f"{counts[g]} foci"))
    sheet(tiles, 5, f"{OUT}/mip_punctate.png")
    json.dump(counts, open(f"{OUT}/spot_counts.json", "w"), indent=1)
    print(counts)


if __name__ == "__main__":
    main()
