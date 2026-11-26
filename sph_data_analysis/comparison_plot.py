"""
Overlay neutrino diagnostics from two (or more) simulations on the same axes.

Reuses the .txt outputs already produced by analyze_neutrinos.py in each model's
analysis directory (no data regeneration).  Both models are drawn on one figure:

  * each quantity (neutrino species / total) keeps the SAME color across models;
  * the models are distinguished by line style (1st = solid, 2nd = dashed, ...);
  * time is measured relative to each model's own merger (single merger marker
    at t - t_merger = 0);
  * the x-range spans all models, each drawn to its own end;
  * a model legend (line style) is placed at the lower right.

With --ejecta-out / --ejecta-time-out it also overlays the models' ejecta
figures fig_ejecta.pdf / fig_ejecta_time.pdf of plot_analysis.py (from the
analyze_ejecta.py outputs), with one color per model and one line style per
ejecta criterion (geodesic solid, Bernoulli dashed).  In the distributions
figure all models are summed up to the last post-merger time they have in
common, re-summed from each model's analysis cache (if a cache is missing,
each model keeps its own window, T_POST_MS in run_analysis.sh).

Example
-------
python comparison_plot.py \
    --dir /lus/.../bundle_run/LR/analysis             --label "eq (corr on)" \
    --dir /lus/.../bundle_run/VVLR_eq_number/analysis --label "non-eq" \
    --out fig_neutrino_comparison.pdf
"""

import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_style  # central house style (repo root); see plot_style.py
from plot_analysis import (load_ejecta_summary, ejecta_cache, ejecta_until,
                           plot_ejecta, plot_ejecta_time)

plot_style.use()

# Unit conversions (match plot_analysis.py)
G          = 6.67430e-11
M_SUN      = 1.98892e30
C          = 2.99792458e8
MSUN_TO_MS = G * M_SUN / C**3 * 1e3   # M_sun -> ms

SPECIES   = ['nue', 'nua', 'nux', 'anux']
SP_LABEL  = [r'$\nu_e$', r'$\bar\nu_e$', r'$\nu_x$', r'$\bar\nu_x$']
SP_COLOR  = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
LINESTYLES = ['-', '--', ':', '-.']
MODEL_COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']   # ejecta figures


def read_merger_time(analysis_dir):
    """Parse t_merger [M_sun] from <analysis_dir>/merger_time.txt (first numeric line)."""
    path = os.path.join(analysis_dir, 'merger_time.txt')
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            return float(s.split()[0])
    raise ValueError(f"No numeric merger time found in {path}")


def load_model(analysis_dir, t_merger=None):
    """Load luminosity + mean-energy series for one model, time relative to merger [ms]."""
    if t_merger is None:
        t_merger = read_merger_time(analysis_dir)

    def load(fname):
        d = np.loadtxt(os.path.join(analysis_dir, fname))
        d[:, 0] = (d[:, 0] - t_merger) * MSUN_TO_MS     # -> t - t_merger [ms]
        return d

    return {
        'Lnu':  {sp: load(f'Lnu_E_{sp}.txt') for sp in SPECIES},
        'Ltot': load('Lnu_E_total.txt'),
        'Eav':  {sp: load(f'Eav_{sp}.txt') for sp in SPECIES},
        't_merger': t_merger,
    }


def span_xlim(models):
    """x-range [ms] about merger covering all models, each to its own end."""
    lefts, rights = [], []
    for m in models:
        t = m['Ltot'][:, 0]
        lefts.append(t.min())
        rights.append(t.max())
    return min(lefts), max(rights)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", action="append", required=True, dest="dirs",
                        help="Analysis directory of a model (repeat for each). "
                             "First is solid, second dashed, etc.")
    parser.add_argument("--label", action="append", dest="labels",
                        help="Legend label for each --dir (repeat, same order).")
    parser.add_argument("--t-merger", action="append", type=float, dest="t_mergers",
                        help="Override merger time [M_sun] per model (else read "
                             "from each model's merger_time.txt).")
    parser.add_argument("--out", default="fig_neutrino_comparison.pdf",
                        help="Output figure filename (written into --outdir).")
    parser.add_argument("--outdir", default=None,
                        help="Where to write the figure "
                             "(default: <first --dir>/comparison_plot/).")
    parser.add_argument("--xlim-ms", nargs=2, type=float, default=None,
                        metavar=("LO", "HI"),
                        help="Override x-range [ms] (default: all models).")
    parser.add_argument("--ejecta-out", default=None,
                        help="Also write the ejecta distributions comparison "
                             "(as fig_ejecta.pdf) to this filename in --outdir.")
    parser.add_argument("--ejecta-time-out", default=None,
                        help="Also write the ejecta history comparison "
                             "(as fig_ejecta_time.pdf) to this filename in --outdir.")
    args = parser.parse_args()

    dirs = args.dirs
    labels = args.labels or [os.path.basename(os.path.normpath(os.path.dirname(d) or d))
                             for d in dirs]
    if len(labels) != len(dirs):
        parser.error("Number of --label must match number of --dir.")
    t_mergers = args.t_mergers if args.t_mergers else [None] * len(dirs)
    if len(t_mergers) != len(dirs):
        parser.error("Number of --t-merger must match number of --dir.")
    if len(dirs) > len(LINESTYLES):
        parser.error(f"At most {len(LINESTYLES)} models supported.")

    models = [load_model(d, tm) for d, tm in zip(dirs, t_mergers)]
    styles = LINESTYLES[:len(models)]

    if args.xlim_ms is not None:
        xlo, xhi = args.xlim_ms
    else:
        xlo, xhi = span_xlim(models)

    fig, axes = plot_style.subplots(1, 2, panel_w=5.5, panel_h=4.2)
    fig.suptitle("BNS — Neutrino emission: model comparison", fontsize=13)

    # ---- Panel 0: energy luminosity (log) ----
    ax = axes[0]
    for m, ls in zip(models, styles):
        ax.semilogy(m['Ltot'][:, 0], m['Ltot'][:, 1], color='black', lw=1.5, ls=ls)
        for sp, col in zip(SPECIES, SP_COLOR):
            d = m['Lnu'][sp]
            ax.semilogy(d[:, 0], d[:, 1], color=col, lw=1.0, ls=ls)
    ax.axvline(0.0, color='gray', ls=':', lw=0.9)
    ax.set_xlabel(r'$t - t_\mathrm{merger}$ [ms]')
    ax.set_ylabel(r'$L_{\nu,E}$  [erg s$^{-1}$]')
    ax.set_ylim(1e51, 3e54)
    ax.set_xlim(xlo, xhi)
    ax.set_title('Energy luminosity')

    # ---- Panel 1: mean neutrino energy (linear) ----
    ax = axes[1]
    for m, ls in zip(models, styles):
        for sp, col in zip(SPECIES, SP_COLOR):
            d = m['Eav'][sp]
            ax.plot(d[:, 0], d[:, 1], color=col, lw=1.2, ls=ls)
    ax.axvline(0.0, color='gray', ls=':', lw=0.9)
    ax.set_xlabel(r'$t - t_\mathrm{merger}$ [ms]')
    ax.set_ylabel(r'$\langle\varepsilon_\nu\rangle$  [MeV]')
    ax.set_ylim(0, 60)
    ax.set_xlim(xlo, xhi)
    ax.set_title('Mean neutrino energy')

    # ---- Legends ----
    # Quantity (color) legend on the luminosity panel.
    qhandles = [Line2D([0], [0], color='black', lw=1.5, label='total')]
    qhandles += [Line2D([0], [0], color=c, lw=1.2, label=l)
                 for c, l in zip(SP_COLOR, SP_LABEL)]
    leg_q = axes[0].legend(handles=qhandles, fontsize=9, ncol=2, loc='upper right')
    axes[0].add_artist(leg_q)

    # Model (line-style) legend at the lower right of each panel.
    mhandles = [Line2D([0], [0], color='dimgray', lw=1.5, ls=ls, label=lab)
                for ls, lab in zip(styles, labels)]
    for ax in axes:
        ax.legend(handles=mhandles, fontsize=9, loc='lower right',
                  framealpha=0.9, title=None)

    fig.tight_layout()
    outdir = args.outdir or os.path.join(dirs[0], 'comparison_plot')
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, args.out)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Saved {out}")
    print(f"  models: " + ", ".join(f"{l} (t_merger={m['t_merger']:g} M_sun, ls='{s}')"
                                     for l, m, s in zip(labels, models, styles)))
    print(f"  x-range: [{xlo:.2f}, {xhi:.2f}] ms about merger")

    if args.ejecta_out or args.ejecta_time_out:
        ejecta = [(lab, col, load_ejecta_summary(d, m['t_merger']))
                  for lab, col, d, m in zip(labels, MODEL_COLORS, dirs, models)]
    if args.ejecta_out:
        # Distributions of all models over the same post-merger window: the
        # last time they all have histograms for, re-summed from their caches.
        caches = [ejecta_cache(d, e, m['t_merger'])
                  for d, (_, _, e), m in zip(dirs, ejecta, models)]
        dist = ejecta
        if all(c is not None for c in caches):
            t_common = min(c['t_end'] for c in caches)
            dist = [(lab, col, ejecta_until(e, c, t_common))
                    for (lab, col, e), c in zip(ejecta, caches)]
            print(f"  ejecta distributions up to the common t - t_merger = {t_common:.2f} ms")
        else:
            missing = [lab for (lab, _, _), c in zip(ejecta, caches) if c is None]
            print("  WARNING: no usable ejecta cache for " + ", ".join(missing)
                  + "; distributions keep each model's own window")
        plot_ejecta(dist, os.path.join(outdir, args.ejecta_out),
                    'BNS SPH — Ejecta mass distributions: model comparison')
    if args.ejecta_time_out:
        plot_ejecta_time(ejecta, os.path.join(outdir, args.ejecta_time_out),
                         'BNS SPH — Ejecta histories: model comparison')


if __name__ == "__main__":
    main()
