"""Simulated fluorescence micrographs of the packed cell (IMPROVEMENTS I2).

Every placed instance that contains the protein carries one fluorophore (a GFP fusion; one per instance, so a homo-
oligomer or complex counts once). Fluorophore positions (instance centres) are binned on a 3D grid, convolved with a
microscope point-spread function, sampled at the camera pixel size and given shot noise, camera offset/read noise and
an autofluorescence background.

PSF: scalar Debye integral for a high-NA objective (aplanatic, sqrt(cos) apodisation; e.g. Gu 2000, Advanced
Optical Imaging Theory, or Richards & Wolf 1959 without polarisation), evaluated per defocus plane:
  h(r, z) = | int_0^alpha sqrt(cos t) J0(k n r sin t) exp(i k n z cos t) sin t dt |^2,  k = 2 pi / lambda, sin alpha = NA/n
Widefield ('wf'): emission PSF (510 nm), all out-of-focus light included. Confocal ('confocal'): excitation (488 nm)
x emission PSF (infinitely small pinhole: an upper bound on sectioning). Defaults: NA 1.4 oil (n 1.515), camera
pixel 6.5 um at 100x = 65 nm. The cell is embedded in index-matched medium (no spherical aberration).

  python microscopy.py GENE [GENE ...] [--z 0] [--out DIR]      -> PNG and .npy of the focal plane
"""
import os, sys, json, argparse
import numpy as np
from scipy import ndimage

from config import WORK

INST = f"{WORK}/contacts/instances.npz"
VOX = 65.0 / 3                # nm, rendering grid (3 x 3 voxels = one 65 nm camera pixel)


class Cell:
    def __init__(self, inst=INST, data=None):
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "simgpu"))
        import boxes
        from contacts import proteins_of
        from config import WEB_OUT
        Z = np.load(inst)
        self.type, self.pos, self.region = Z["type"].astype(np.int64), Z["pos"].astype(np.float64), Z["region"]
        self.types = boxes.Library(data or WEB_OUT).types
        names, off, mem = proteins_of(self.types)
        self.proteins = names
        self.by_protein = {}
        for t in range(len(self.types)):
            for k in mem[off[t]:off[t + 1]]:
                self.by_protein.setdefault(names[k], []).append(t)
        self.lo = self.pos.min(0) - 400; self.hi = self.pos.max(0) + 400
        self.shape = tuple(np.ceil((self.hi - self.lo) / VOX).astype(int))

    def positions(self, gene):
        ts = self.by_protein.get(gene)
        if not ts:
            raise KeyError(gene)
        m = np.isin(self.type, ts)
        return self.pos[m]

    def density(self, pts):
        idx = np.floor((pts - self.lo) / VOX).astype(int)
        g = np.zeros(self.shape, np.float32)
        np.add.at(g, tuple(idx.T), 1.0)
        return g


def debye_psf(lam, na=1.4, n=1.515, rmax=2000.0, zs=(0.0,), dr=VOX):
    """radial PSF profiles h(r) for each defocus z (nm); normalised so that the in-focus plane integrates to 1"""
    from scipy.special import j0
    alpha = np.arcsin(na / n)
    t = np.linspace(0, alpha, 400)
    r = np.arange(0, rmax + dr, dr / 2)
    k = 2 * np.pi / lam
    J = j0(k * n * np.outer(r, np.sin(t)))                       # (nr, nt)
    w = np.sqrt(np.cos(t)) * np.sin(t)
    out = []
    for z in zs:
        f = J * (w * np.exp(1j * k * n * z * np.cos(t)))[None, :]
        out.append(np.abs(np.trapezoid(f, t, axis=1)) ** 2)
    out = np.array(out)
    norm = np.trapezoid(out[np.argmin(np.abs(np.asarray(zs)))] * 2 * np.pi * r, r)
    return r, out / norm


def kernel2d(r, h, half):
    """2D kernel on the VOX grid from a radial profile"""
    x = np.arange(-half, half + 1) * VOX
    R = np.hypot(*np.meshgrid(x, x, indexing="ij"))
    return np.interp(R, r, h, right=0.0) * VOX * VOX


def render(cell, pts, axis=2, z=0.0, na=1.4, pixel=65.0, mode="wf", zrange=2000.0, photons=None, n_imm=1.515,
           background=2.0, offset=100.0, read_noise=1.6, seed=0, rmax=1500.0):
    """focal plane image normal to `axis` at coordinate z (nm). photons: mean detected photons per fluorophore in
    focus (None -> noiseless, in fluorophore units)"""
    from scipy.signal import fftconvolve
    g = cell.density(pts)
    k0 = int(np.floor((z - cell.lo[axis]) / VOX))
    zr = int(np.ceil(zrange / VOX))
    sl = [slice(None)] * 3; sl[axis] = slice(max(0, k0 - zr), min(g.shape[axis], k0 + zr + 1))
    sub = np.moveaxis(g[tuple(sl)], axis, 0)
    dz = (np.arange(sub.shape[0]) + max(0, k0 - zr) - k0) * VOX
    r, he = debye_psf(510.0, na, n=n_imm, zs=dz, rmax=rmax)
    if mode == "confocal":
        _, hx = debye_psf(488.0, na, n=n_imm, zs=dz, rmax=rmax)
        h = he * hx
        h /= np.trapezoid(h[np.argmin(np.abs(dz))] * 2 * np.pi * r, r)
    else:
        h = he
    half = int(rmax / VOX)
    img = np.zeros(sub.shape[1:], np.float64)
    for i in range(sub.shape[0]):
        if sub[i].any():
            img += fftconvolve(sub[i].astype(np.float64), kernel2d(r, h[i], half), mode="same")
    f = max(1, int(round(pixel / VOX)))
    hh, ww = (img.shape[0] // f) * f, (img.shape[1] // f) * f
    img = img[:hh, :ww].reshape(hh // f, f, ww // f, f).sum((1, 3))
    if photons is None:
        return img
    rng = np.random.default_rng(seed)
    lam_px = np.clip(img, 0, None) * photons + background
    return rng.poisson(lam_px) + offset + rng.normal(0, read_noise, lam_px.shape)


def render_stack(cell, pts, step=200.0, axis=2, **kw):
    """z-stack through the whole cell (planes every `step` nm): array (nz, ny, nx) and plane coordinates"""
    z = np.arange(cell.pos[:, axis].min() + step / 2, cell.pos[:, axis].max(), step)
    return np.array([render(cell, pts, axis=axis, z=zz, **kw) for zz in z]), z


def count_spots(stack, cellmask=None, sigma_px=1.0, contrast=0.3, min_sep_px=3, ring=(4, 7)):
    """foci in a z-stack: local maxima of the (Gaussian-smoothed) max projection inside the cell that are at least
    `contrast` (fractional) brighter than the median of a surrounding annulus (radii `ring` px, 65 nm px) -- a
    background-relative criterion, so diffuse signal with noise texture gives few spots. Returns (n, coordinates)."""
    mip = ndimage.gaussian_filter(stack.max(0).astype(float), sigma_px)
    mx = ndimage.maximum_filter(mip, size=2 * min_sep_px + 1)
    inside = np.ones_like(mip, bool) if cellmask is None else cellmask
    cand = np.argwhere((mip == mx) & inside)
    yy, xx = np.mgrid[-ring[1]:ring[1] + 1, -ring[1]:ring[1] + 1]
    rr = np.hypot(yy, xx); ann = (rr >= ring[0]) & (rr <= ring[1])
    pad = np.pad(mip, ring[1], mode="edge")
    keep = []
    for (i, j) in cand:
        bg = np.median(pad[i:i + 2 * ring[1] + 1, j:j + 2 * ring[1] + 1][ann])
        if mip[i, j] > (1 + contrast) * bg and mip[i, j] - bg > 0.05 * (mip[inside].max() - np.median(mip[inside])):
            keep.append((i, j))
    return len(keep), np.array(keep)


def to_png(img, path, lo_pct=1, hi_pct=99.7):
    from PIL import Image
    a = np.asarray(img, float)
    lo, hi = np.percentile(a, lo_pct), np.percentile(a, hi_pct)
    u8 = (np.clip((a - lo) / max(hi - lo, 1e-9), 0, 1) ** 0.8 * 255).astype(np.uint8)
    Image.fromarray(u8.T[::-1]).save(path)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("genes", nargs="+"); ap.add_argument("--z", type=float, default=0.0); ap.add_argument("--axis", type=int, default=2)
    ap.add_argument("--out", default=f"{WORK}/microscopy"); ap.add_argument("--photons", type=float, default=None)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    c = Cell()
    for gname in a.genes:
        im = render(c, c.positions(gname), axis=a.axis, z=a.z, photons=a.photons)
        np.save(f"{a.out}/{gname}_z{int(a.z)}.npy", im); to_png(im, f"{a.out}/{gname}_z{int(a.z)}.png")
        print(gname, im.shape, f"{im.sum():.0f}")
