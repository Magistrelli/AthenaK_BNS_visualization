"""
Opacity / radiation slice COMPARISON between two AthenaK runs.

For each selected variable and each matched output time, renders a 3-panel
figure:  [run A]  [run B]  [ratio A/B] , on the xz slice.  The run panels use
the same colormap/norm/units as BNS_all.py (via athenak_units.plot_params);
the ratio panel is a diverging log map centered on 1, so you can see at a glance
*where* and *by how much* the two runs' opacities differ.

The two runs' AMR grids may differ, so each is rasterized onto a common uniform
grid (nearest-cell) before the ratio is taken.

Conventions mirror BNS_all.py (xz slice, M_sun or km axes, AMR-aware) and the
two-folder comparison style of sph_data_analysis/run_comparison.sh.

Example
-------
python BNS_compare.py \
    --run /lus/.../VVLR_test_2 --label "eq" \
    --run /lus/.../VVLR_test_3 --label "non-eq" \
    --vars eta_1:2 abs_1:2 eta_0:2 eta_1:0 \
    --units cgs --outdir /lus/.../VVLR_test_2/compare_xz
"""

import argparse
import glob
import os
import sys
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, SymLogNorm

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.environ.get("ATHPLOT_PATH", os.path.join(_HERE, "plot-tools")))
sys.path.append(_HERE)
from athplot.load_ath_bin import BinaryData   # noqa: E402
import athenak_units as au                      # noqa: E402
import plot_style                               # noqa: E402

plot_style.use()  # central house style (see plot_style.py)

from BNS_all import combine_hst

# variable shorthand -> per-plane (new dump prefix, old prefix)
# NB: the "old" fallback (rad_m1_*) is the SPHERICAL-surface dump, not a slice,
# so it is only a valid fallback for xz-era runs.  For xy we require a genuine
# *_xy slice dump (no spherical fallback), else the plane silently reports none.
VAR_PREFIX = {
    'eta_0': {'xz': ("bns.m1opac_xz", "bns.m1_opac_xz", "bns.eta_0_xz", "bns.rad_m1_eta_0"),
              'xy': ("bns.m1opac_xy", "bns.m1_opac_xy", "bns.eta_0_xy")},
    'eta_1': {'xz': ("bns.m1opac_xz", "bns.m1_opac_xz", "bns.eta_1_xz", "bns.rad_m1_eta_1"),
              'xy': ("bns.m1opac_xy", "bns.m1_opac_xy", "bns.eta_1_xy")},
    'abs_0': {'xz': ("bns.m1opac_xz", "bns.m1_opac_xz", "bns.abs_0_xz", "bns.rad_m1_abs_0"),
              'xy': ("bns.m1opac_xy", "bns.m1_opac_xy", "bns.abs_0_xy")},
    'abs_1': {'xz': ("bns.m1opac_xz", "bns.m1_opac_xz", "bns.abs_1_xz", "bns.rad_m1_abs_1"),
              'xy': ("bns.m1opac_xy", "bns.m1_opac_xy", "bns.abs_1_xy")},
    'E':     {'xz': ("bns.m1_xz", "bns.m1_mom_xz", "bns.rad_E_xz", "bns.rad_m1_E"),
              'xy': ("bns.m1_xy", "bns.m1_mom_xy", "bns.rad_E_xy")},
    'N':     {'xz': ("bns.m1_xz", "bns.m1_mom_xz", "bns.rad_N_xz", "bns.rad_m1_N"),
              'xy': ("bns.m1_xy", "bns.m1_mom_xy", "bns.rad_N_xy")},
    # fluid fields live in the prim dump (no species suffix)
    'temperature': {'xz': ("bns.prim_xz", "bns.prim_xz"),
                    'xy': ("bns.prim_xy", "bns.prim_xy")},
    's_00':        {'xz': ("bns.prim_xz", "bns.prim_xz"),
                    'xy': ("bns.prim_xy", "bns.prim_xy")},
}
SP_LABEL = {0: r"\nu_e", 1: r"\bar\nu_e", 2: r"\nu_x", 3: r"\bar\nu_x"}
# colormap for the run panels — match BNS_all.py style (rates=inferno, T=viridis, Ye=plasma)
RATE_CMAP = "inferno"
VAR_CMAP = {'temperature': 'viridis', 's_00': 'plasma'}
# fixed display range (linear), matching BNS_all.py; overrides auto-scaling
VAR_RANGE = {'temperature': (0.0, 70.0), 's_00': (0.0, 0.5)}
# third panel = SIGNED difference (A-B) on a symlog scale spanning +/- this
# magnitude (diverging, centered on 0):
VAR_DIFF = {'temperature': 3.0, 's_00': 0.1}
DIFF_UNIT = {'temperature': r"\mathrm{MeV}", 's_00': r""}


def spatial_units(units):
    """Axis scale (M_sun -> displayed) and label; mirrors BNS_all.spatial_units."""
    if units in ("cgs", "ngs"):
        return 1.476625038, r"\mathrm{km}"   # 1 M_sun in km (G=c=1)
    return 1.0, r"\mathrm{M}_\odot"


def find_output_dirs(simpath):
    ds = sorted(glob.glob(os.path.join(simpath, "output-*")))
    return ds if ds else [simpath]


def _normalise_prefixes(prefixes):
    if isinstance(prefixes, str):
        return (prefixes,)
    return tuple(prefixes)


def detect_prefix(output_dirs, *prefixes):
    """Keep candidates ordered; frames_by_time resolves each segment."""
    if not prefixes:
        raise ValueError("at least one output prefix is required")
    return tuple(prefixes)


_FRAME_CACHE = {}


def frames_by_time(output_dirs, prefixes):
    """{rounded_time: filepath}; newest naming is selected per segment."""
    prefixes = _normalise_prefixes(prefixes)
    key = (prefixes, tuple(output_dirs))
    if key in _FRAME_CACHE:
        return _FRAME_CACHE[key]
    tm = {}
    for d in output_dirs:
        for prefix in prefixes:
            files = sorted(glob.glob(os.path.join(d, "bin", prefix + ".*.bin")))
            if not files:
                continue
            for f in files:
                try:
                    t = BinaryData(f).time
                except Exception:
                    continue
                if np.isfinite(t):
                    tm[round(t, 1)] = f
            break
    _FRAME_CACHE[key] = tm
    return tm


def rasterize(bd, var, xlim, zlim, nx, nz):
    """Sample an AthenaK xz slice onto a uniform (nz, nx) grid (nearest cell)."""
    xs = np.linspace(xlim[0], xlim[1], nx)
    zs = np.linspace(zlim[0], zlim[1], nz)
    grid = np.full((nz, nx), np.nan)
    for block, data in bd.get_block_data(var):
        x0, x1, z0, z1 = block.get_extent()
        nbz, nbx = data.shape
        dxb = (x1 - x0) / nbx
        dzb = (z1 - z0) / nbz
        cx = np.where((xs >= x0) & (xs < x1))[0]
        cz = np.where((zs >= z0) & (zs < z1))[0]
        if cx.size == 0 or cz.size == 0:
            continue
        ibx = np.clip(((xs[cx] - x0) / dxb).astype(int), 0, nbx - 1)
        ibz = np.clip(((zs[cz] - z0) / dzb).astype(int), 0, nbz - 1)
        grid[np.ix_(cz, cx)] = data[np.ix_(ibz, ibx)]
    return xs, zs, grid


def render_one(var, t, fA, fB, labels, units, outdir, out_prefix,
               xlim, zlim, nx, nz, rmax, coord_scale, coord_label, vlabel,
               rate_vmin=None, rate_vmax=None):
    has_sp = ":" in var
    base, _, sp = var.partition(":")
    sp = int(sp) if sp else 0
    pp_key = f"{base}:0" if has_sp else base   # opacity params keyed on species 0
    scale, plabel, use_log, vmin, vmax = au.plot_params(pp_key, units)
    cmap = VAR_CMAP.get(base, RATE_CMAP)
    if base in VAR_RANGE:
        # fluid fields (T, Ye): fixed BNS_all.py display range, linear
        vmin, vmax = VAR_RANGE[base]
    else:
        # rate color range: prefer explicit / auto-tuned range over the very wide
        # plot_params defaults (which wash the rates out to one color)
        if rate_vmin is not None:
            vmin = rate_vmin
        if rate_vmax is not None:
            vmax = rate_vmax

    bdA, bdB = BinaryData(fA), BinaryData(fB)
    _, _, gA = rasterize(bdA, var, xlim, zlim, nx, nz)
    _, _, gB = rasterize(bdB, var, xlim, zlim, nx, nz)
    gA = gA * scale
    gB = gB * scale

    ext = [xlim[0] * coord_scale, xlim[1] * coord_scale,
           zlim[0] * coord_scale, zlim[1] * coord_scale]
    run_norm = (LogNorm(vmin=vmin, vmax=vmax) if use_log
                else Normalize(vmin=vmin, vmax=vmax))

    is_diff = base in VAR_DIFF
    if is_diff:
        # third panel = SIGNED difference (A - B), diverging symlog
        third = gA - gB
        dmax = VAR_DIFF[base]
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            third = np.where((gA > 0) & (gB > 0), gA / gB, np.nan)
            # mask cells where BOTH runs are below the colorbar floor (negligible
            # signal) — otherwise the ratio of two near-zero values is pure noise.
            floor = vmin if (vmin and vmin > 0) else 0.0
            third = np.where((gA > floor) | (gB > floor), third, np.nan)
            # clip into the displayed decade so out-of-range cells saturate
            third = np.clip(third, 1.0 / rmax, rmax)

    fig, axes = plot_style.subplots(1, 3, panel_w=4.5, panel_h=4.6,
                                    constrained_layout=True)
    _NICE = {'temperature': r"T~[\mathrm{MeV}]", 's_00': r"Y_e"}
    if has_sp:
        title = rf"${base}\;[{SP_LABEL.get(sp, str(sp))}]$"
    else:
        title = rf"${_NICE.get(base, base)}$"
    fig.suptitle(rf"{title}   $t={t:.1f}\,\mathrm{{M}}_\odot$", fontsize=13)

    for ax, g, ttl in ((axes[0], gA, labels[0]), (axes[1], gB, labels[1])):
        im = ax.imshow(g, origin="lower", extent=ext, cmap=cmap,
                       norm=run_norm, interpolation="nearest")
        ax.set_title(ttl)
        ax.set_xlabel(rf"$x~[{coord_label}]$")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label=plabel)
    axes[0].set_ylabel(rf"${vlabel}~[{coord_label}]$")

    if is_diff:
        # signed-difference panel: diverging LINEAR, +/- dmax about 0
        du = DIFF_UNIT.get(base, "")
        _dcmap = matplotlib.colormaps["RdBu_r"].copy()
        im = axes[2].imshow(third, origin="lower", extent=ext, cmap=_dcmap,
                            norm=Normalize(vmin=-dmax, vmax=dmax),
                            interpolation="nearest")
        axes[2].set_title(f"{labels[0]} $-$ {labels[1]}")
        axes[2].set_xlabel(rf"$x~[{coord_label}]$")
        ticks = [-dmax, -0.5 * dmax, 0.0, 0.5 * dmax, dmax]
        _dl = rf"$\Delta~[{du}]$" if du else r"$\Delta Y_e$"
        cb = fig.colorbar(im, ax=axes[2], fraction=0.046, pad=0.02, ticks=ticks,
                          label=_dl)
        cb.ax.set_yticklabels([f"{v:g}" for v in ticks])
    else:
        # ratio panel: diverging map, LOG scale centered on 1 (symmetric: 1/rmax..rmax)
        _rcmap = matplotlib.colormaps["RdBu_r"].copy()
        _rcmap.set_bad("0.85")   # negligible-signal cells shown as light gray
        im = axes[2].imshow(third, origin="lower", extent=ext, cmap=_rcmap,
                            norm=LogNorm(vmin=1.0 / rmax, vmax=rmax),
                            interpolation="nearest")
        axes[2].set_title(f"{labels[0]} / {labels[1]}")
        axes[2].set_xlabel(rf"$x~[{coord_label}]$")
        # symmetric ticks in log about 1, e.g. 0.1, 0.316, 1, 3.16, 10 for rmax=10
        ticks = sorted({1.0 / rmax, 1.0 / np.sqrt(rmax), 1.0,
                        np.sqrt(rmax), rmax})
        cb = fig.colorbar(im, ax=axes[2], fraction=0.046, pad=0.02, ticks=ticks,
                          label=f"ratio  (white=equal, blue<1<red, sat @ {rmax:g}$\\times$)")
        cb.ax.set_yticklabels([f"{v:g}" for v in ticks])

    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{out_prefix}_{t:09.1f}.png")
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    med = np.nanmedian(third)
    return out, med


def read_user_header(output_dir):
    """Read column names from bns.user.hst header."""

    outdir = output_dir[0] if isinstance(output_dir, list) else output_dir
    user_file = os.path.join(outdir, "bns.user.hst")
    if not os.path.isfile(user_file): return None
    names = []

    with open(user_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                names += re.findall(r"\[\d+\]=([^\s]+)", line)

    return names


def compare_hst_all(outdirs, labels, figpath):
    """
    Compare history variables from two runs in a single figure.

    Panels:
      1) Baryonic mass (bns.mhd.hst)
      2) Hamiltonian constraint norm (bns.z4c.user.hst)
      3+) All quantities from bns.user.hst (header-driven)
    """

    os.makedirs(figpath, exist_ok=True)

    # Read histories
    mhd_A = combine_hst(outdirs[0], "bns.mhd.hst")
    mhd_B = combine_hst(outdirs[1], "bns.mhd.hst")

    z4c_A = combine_hst(outdirs[0], "bns.z4c.user.hst")
    z4c_B = combine_hst(outdirs[1], "bns.z4c.user.hst")

    user_A = combine_hst(outdirs[0], "bns.user.hst")
    user_B = combine_hst(outdirs[1], "bns.user.hst")

    # Column of each user quantity, per run, from that run's own header
    # (skip time and dt).  A quantity missing in one run is plotted for the
    # other run only.
    user_cols = []
    for outdir, user in zip(outdirs, (user_A, user_B)):
        ncol = user.shape[1] if user is not None else 0
        names = read_user_header(outdir) or [f"user[{i}]" for i in range(ncol)]
        user_cols.append({n: i for i, n in enumerate(names) if 2 <= i < ncol})
    user_quantities = list(user_cols[0]) + [n for n in user_cols[1]
                                            if n not in user_cols[0]]

    def plot_col(ax, data, col, label):
        """Plot column `col` vs time, if this run has it."""
        if data is not None and col < data.shape[1]:
            ax.plot(data[:, 0], data[:, col], label=label)

    # Number of panels:
    # mhd + z4c + all user quantities
    nplots = 2 + len(user_quantities)
    ncols = 2
    nrows = int(np.ceil(nplots / ncols))
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(12, 4*nrows),
        constrained_layout=True
    )
    axes = np.atleast_1d(axes).flatten()

    # Baryonic mass
    ax = axes[0]
    plot_col(ax, mhd_A, 2, labels[0])
    plot_col(ax, mhd_B, 2, labels[1])
    ax.set_xlabel("t (Msun)")
    ax.set_ylabel("Baryonic mass")
    ax.legend()

    # Hamiltonian constraint
    ax = axes[1]
    plot_col(ax, z4c_A, 3, labels[0])
    plot_col(ax, z4c_B, 3, labels[1])
    ax.set_xlabel("t (Msun)")
    ax.set_ylabel("H-norm2")
    ax.set_yscale("log")
    ax.legend()

    # User history: use header names
    for i, name in enumerate(user_quantities):
        ax = axes[i + 2]
        for user, cols, label in zip((user_A, user_B), user_cols, labels):
            if name in cols:
                plot_col(ax, user, cols[name], label)

        ax.set_xlabel("t (Msun)")
        ax.set_ylabel(name)
        ax.legend()

    # Remove unused axes
    for ax in axes[nplots:]:
        ax.remove()

    fig.suptitle("History comparison")

    outfile = os.path.join(figpath, "history_comparison.png")
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"  Saved {outfile}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, dest="runs",
                    help="Run directory (contains output-*/bin/). Give exactly two.")
    ap.add_argument("--label", action="append", dest="labels",
                    help="Legend label per --run (same order).")
    ap.add_argument("--vars", nargs="+",
                    default=["eta_1:2", "abs_1:2", "eta_0:2", "eta_1:0", "abs_1:0"],
                    help="Variables as base:species, e.g. eta_1:2 abs_1:2 E:0.")
    ap.add_argument("--no-hst", action="store_true",
                    help="Do not generate comparison plots of history variables")
    ap.add_argument("--units", default="cgs", choices=["code", "cgs", "ngs"])
    ap.add_argument("--plane", default="xz", choices=["xz", "xy"],
                    help="Slice plane to compare (needs the matching *_<plane> dumps).")
    ap.add_argument("--xlim", nargs=2, type=float, default=[-60, 60],
                    help="x-range [M_sun].")
    ap.add_argument("--zlim", nargs=2, type=float, default=[-60, 60],
                    help="vertical-axis range [M_sun] (z for xz, y for xy).")
    ap.add_argument("--nx", type=int, default=400)
    ap.add_argument("--nz", type=int, default=400)
    ap.add_argument("--rmax", type=float, default=10.0,
                    help="Ratio color half-range: panel spans 1/rmax .. rmax "
                         "(log, symmetric about 1). Default 10 -> 0.1 .. 1 .. 10.")
    ap.add_argument("--rate-vmin", type=float, default=None,
                    help="Override rate (left panels) color min [display units].")
    ap.add_argument("--rate-vmax", type=float, default=None,
                    help="Override rate (left panels) color max [display units].")
    ap.add_argument("--outdir", default=None,
                    help="Output dir (default: <run A>/compare_<plane>/).")
    ap.add_argument("--every", type=int, default=1,
                    help="Use every Nth matched time (thin frames).")
    args = ap.parse_args()

    if len(args.runs) != 2:
        ap.error("Provide exactly two --run.")
    labels = args.labels or [os.path.basename(os.path.normpath(r)) for r in args.runs]
    if len(labels) != 2:
        ap.error("Need one --label per --run.")
    coord_scale, coord_label = spatial_units(args.units)
    vlabel = "y" if args.plane == "xy" else "z"
    outdir = args.outdir or os.path.join(args.runs[0],
                                         f"compare_{args.plane}")
    odirs = [find_output_dirs(r) for r in args.runs]

    print(f"A = {labels[0]} : {args.runs[0]}")
    print(f"B = {labels[1]} : {args.runs[1]}")
    print(f"plane={args.plane}  units={args.units}  xlim={args.xlim}  "
          f"zlim={args.zlim}  rmax={args.rmax}")
    print(f"outdir = {outdir}\n")

    if not args.no_hst:
        compare_hst_all(odirs, labels, outdir)

    for var in args.vars:
        base = var.partition(":")[0]
        if base not in VAR_PREFIX:
            print(f"[skip] unknown var base '{base}' in '{var}'"); continue
        prefA = detect_prefix(odirs[0], *VAR_PREFIX[base][args.plane])
        prefB = detect_prefix(odirs[1], *VAR_PREFIX[base][args.plane])
        tA = frames_by_time(odirs[0], prefA)
        tB = frames_by_time(odirs[1], prefB)
        common = sorted(set(tA) & set(tB))[::args.every]
        if not common:
            print(f"[{var}] no common times (A:{len(tA)} B:{len(tB)}) "
                  f"for plane={args.plane} (prefix {prefA})"); continue
        out_prefix = f"cmp_{var.replace(':', '_')}"
        # auto rate range from the LAST (most active) matched frame unless
        # explicitly overridden -> fixes the "one color everywhere" washout
        rv_min, rv_max = args.rate_vmin, args.rate_vmax
        if base in VAR_RANGE:
            rv_min, rv_max = VAR_RANGE[base]     # fixed BNS_all.py range
        elif rv_min is None or rv_max is None:
            scale, _, use_log, _, _ = au.plot_params(f"{base}:0", args.units)
            g = []
            for f in (tA[common[-1]], tB[common[-1]]):
                _, _, gg = rasterize(BinaryData(f), var, args.xlim, args.zlim,
                                     args.nx, args.nz)
                g.append((gg * scale).ravel())
            g = np.concatenate(g)
            g = g[np.isfinite(g) & (g > 0)]
            if g.size:
                if use_log:
                    # top ~5 decades below the peak -> good contrast in the
                    # bright (emitting) region; faint outskirts saturate low
                    hi = 10 ** np.ceil(np.log10(np.percentile(g, 99.8)))
                    if rv_max is None:
                        rv_max = hi
                    if rv_min is None:
                        rv_min = rv_max / 1e5
                else:
                    lo, hi = np.percentile(g, [2.0, 99.8])
                    if rv_min is None:
                        rv_min = lo
                    if rv_max is None:
                        rv_max = hi
            else:
                # Quantity is zero everywhere (LogNorm needs vmin > 0).
                if rv_max is None:
                    rv_max = 1.0
                if rv_min is None:
                    rv_min = rv_max / 1e5 if use_log else 0.0

        print(f"[{var}] {len(common)} matched times "
              f"[{common[0]:.0f},{common[-1]:.0f}] M_sun  rate range=[{rv_min:.2e},{rv_max:.2e}]")
        for i, t in enumerate(common):
            try:
                _, med = render_one(var, t, tA[t], tB[t], labels, args.units,
                                    outdir, out_prefix, args.xlim, args.zlim,
                                    args.nx, args.nz, args.rmax,
                                    coord_scale, coord_label, vlabel,
                                    rate_vmin=rv_min, rate_vmax=rv_max)
            except Exception as e:
                print(f"  t={t:.1f} FAILED: {e}"); continue
            if i == 0 or i == len(common) - 1:
                print(f"  t={t:.1f}  median(A/B)={med:.3f}")
        print(f"  -> {outdir}/{out_prefix}_*.png")


if __name__ == "__main__":
    main()
