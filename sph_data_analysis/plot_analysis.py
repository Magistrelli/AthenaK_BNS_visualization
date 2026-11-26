"""
Plot ejecta and neutrino diagnostics produced by analyze_ejecta.py and
analyze_neutrinos.py.

Reads the .txt outputs from --output-dir and writes:
    fig_ejecta.pdf          (1x5 mass distributions: v_inf, T, Y_e, entropy, angle)
    fig_ejecta_time.pdf     (1x3 ejecta time evolution: rate, mass, <Y_e>)
    fig_ejecta_ye_theta.pdf (2-D Y_e vs theta composition map)
    fig_neutrino.pdf        (2 panels: neutrino luminosity + mean energy)
"""

import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_style  # central house style (repo root); see plot_style.py

plot_style.use()

# Unit conversions
G       = 6.67430e-11
M_SUN   = 1.98892e30
C       = 2.99792458e8
MSUN_TO_MS = G * M_SUN / C**3 * 1e3


def cumulative(arr):
    t = arr[:, 0]; r = arr[:, 1]; dt = np.diff(t)
    return t[1:], np.cumsum(0.5 * (r[:-1] + r[1:]) * dt)


# Ejecta panels with a log y-axis show this many decades below their peak;
# the near-empty pre-merger / high-v tails would otherwise set the range.
LOG_DECADES = 4


def _logy(ax, *ys):
    """Log y-axis spanning LOG_DECADES below the largest positive value."""
    ax.set_yscale('log')
    top = max((np.max(y[y > 0]) for y in ys if np.any(y > 0)), default=None)
    if top is not None:
        ax.set_ylim(top * 10**-LOG_DECADES, top * 2)


def _fmt_mass(m):
    return f'{m:.2e}'.replace('e-0', 'e-').replace('e+0', 'e')


# Ejecta criteria: one line style each.  A single model also gets one color per
# criterion; in a model comparison the color tells the models apart.
CRIT_NAME  = {'geo': 'geodesic', 'Ber': 'Bernoulli'}
CRIT_COLOR = {'geo': 'tab:blue', 'Ber': 'tab:orange'}
CRIT_LS    = {'geo': '-', 'Ber': '--'}


def _hist_window(d):
    """t - t_merger [ms] up to which the distributions of `d` were summed.

    The window (T_POST_MS in run_analysis.sh) is not stored in the outputs: it
    is where the time-integrated ejection rate reaches the distributions'
    total mass, or the end of the run if it never does.
    """
    t, M = cumulative(d['rate'])
    i = np.searchsorted(M, d['ye'][:, 1].sum() * (1 - 1e-3))
    return t[min(i, len(t) - 1)]


def _window_text(models):
    """Post-merger window of the distributions, for the figure title."""
    ends = [(label, m['t_hist']) for label, _, m in models]
    text = r'distributions up to $t - t_\mathrm{merger}$ = '
    if len({round(t) for _, t in ends}) == 1:
        return text + f'{ends[0][1]:.0f} ms'
    return text + ', '.join(f'{t:.0f} ms ({label})' for label, t in ends)


def load_ejecta_summary(adir, t_merger):
    """Inputs of the ejecta figures, per criterion ('geo', 'Ber').

    Time series are (t - t_merger [ms], value), with the ejection rate in
    M_sun/ms; the distributions cover the post-merger histogram window.
    """
    def load(f): return np.loadtxt(os.path.join(adir, f))

    def t_rel(a):
        a = a.copy()
        a[:, 0] = (a[:, 0] - t_merger) * MSUN_TO_MS
        return a

    out = {}
    for c in CRIT_NAME:
        rate = t_rel(load(f'Mej_rate_{c}.txt'))
        rate[:, 1] /= MSUN_TO_MS
        out[c] = {'rate':   rate,
                  'ye_avg': t_rel(load(f'Ye_avg_{c}.txt')),
                  'vinf':   load(f'Mej_vinf_{c}.txt'),
                  'T':      load(f'Mej_T_{c}.txt'),
                  'ye':     load(f'Mej_Ye_{c}.txt'),
                  's':      load(f'Mej_s_{c}.txt'),
                  'rhoej':  load(f'rhoej_theta_{c}.txt')}
    out['t_hist'] = _hist_window(out['geo'])
    return out


def _ejecta_lines(models):
    """(data, color, line style, legend name) of every curve of an ejecta figure.

    ``models`` is a list of (label, color, load_ejecta_summary(...)); a single
    model is passed as (None, None, summary).
    """
    multi = len(models) > 1
    for label, mcol, m in models:
        for c in CRIT_NAME:
            name = f'{label} {CRIT_NAME[c]}' if multi else CRIT_NAME[c]
            yield m[c], mcol or CRIT_COLOR[c], CRIT_LS[c], name


def plot_ejecta(models, out, title):
    """Ejecta mass distributions vs v_inf, temperature, Y_e, entropy and angle.

    One row of five panels; the first four show ejected mass per bin [M_sun]
    on a log scale fixed to 1e-6 .. 1e-2, the last dM/dtheta.  The geodesic
    and Bernoulli criteria are overlaid.  ``models`` as for _ejecta_lines.
    """
    # (summary key, x-label, fixed xlim or None -> auto to populated bins,
    #  panel title, optional vertical guide line)
    panels = [
        ('vinf', r'$v_\infty$ [$c$]', (0, 0.8), 'Velocity', None),
        ('T', r'$T$ [MeV]', None, 'Temperature', None),
        ('ye', r'$Y_e$', (0, 0.6), r'$Y_e$', 0.5),
        ('s', r'$s$ [$k_B\,\mathrm{baryon}^{-1}$]', None, 'Entropy', None),
    ]

    fig, axes = plot_style.subplots(1, 5, panel_w=4.2, panel_h=4.0)
    fig.suptitle(f'{title}\n{_window_text(models)}', fontsize=13)

    for ax, (key, xlabel, xlim, ptitle, vline) in zip(axes, panels):
        xmax = 0.0
        for d, col, ls, name in _ejecta_lines(models):
            h = d[key]
            ax.step(h[:, 0], h[:, 1], color=col, ls=ls, where='mid',
                    lw=1.2, label=name)
            # auto-focus on the populated bins (all curves)
            nz = np.where(h[:, 1] > 0)[0]
            xmax = max(xmax, h[nz[-1], 0] * 1.1 if len(nz) else h[-1, 0])
        if vline is not None:
            ax.axvline(vline, color='gray', ls=':', lw=0.8)
        ax.set_yscale('log')
        ax.set_ylim(1e-6, 1e-2)
        ax.set_xlabel(xlabel)
        ax.set_title(ptitle)
        ax.set_xlim(*(xlim if xlim is not None else (0, xmax)))

    ax = axes[4]
    ys = []
    for d, col, ls, name in _ejecta_lines(models):
        h = d['rhoej']
        ax.plot(h[:, 0] * 180 / np.pi, h[:, 1], color=col, ls=ls, label=name, lw=1.2)
        ys.append(h[:, 1])
    _logy(ax, *ys)
    ax.set_xlabel(r'$\theta$ [deg]')
    ax.set_ylabel(r'$dM_\mathrm{ej}/d\theta$  [$M_\odot\,\mathrm{rad}^{-1}$]')
    ax.set_xlim(0, 180); ax.xaxis.set_major_locator(MultipleLocator(30))
    ax.set_title('Angle')

    axes[0].set_ylabel(r'$M_\mathrm{ej}$ per bin  [$M_\odot$]')
    axes[0].legend(fontsize=8 if len(models) > 1 else 9)

    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved {out}')


def plot_ejecta_time(models, out, title):
    """Ejecta histories: ejection rate, cumulative mass and <Ye>.

    Masses in M_sun; times relative to the merger.  ``models`` as for
    _ejecta_lines.  Every model's time series runs to its own end, and its
    quoted totals are taken there.
    """
    multi = len(models) > 1
    xlab = r'$t - t_\mathrm{merger}$ [ms]'
    xlim_t = dict(left=-2.0)
    legend_fs = 8 if multi else 9

    fig, axes = plot_style.subplots(1, 3, panel_w=5.0, panel_h=4.0)
    fig.suptitle(title, fontsize=13)

    ax = axes[0]
    ys = []
    for d, col, ls, name in _ejecta_lines(models):
        r = d['rate']
        ax.plot(r[:, 0], r[:, 1], color=col, ls=ls, label=name, lw=1.2)
        ys.append(r[:, 1])
    _logy(ax, *ys)
    ax.set_xlabel(xlab); ax.set_ylabel(r'$\dot{M}_\mathrm{ej}$  [$M_\odot\,\mathrm{ms}^{-1}$]')
    ax.legend(fontsize=legend_fs); ax.set_xlim(**xlim_t); ax.set_title('Mass ejection rate')

    ax = axes[1]
    ys = []
    for d, col, ls, name in _ejecta_lines(models):
        t, M = cumulative(d['rate'])
        ax.plot(t, M, color=col, ls=ls, lw=1.2,
                label=name if multi else f'{name} ({_fmt_mass(M[-1])} $M_\\odot$)')
        ys.append(M)
    _logy(ax, *ys)
    ax.set_xlabel(xlab); ax.set_ylabel(r'$M_\mathrm{ej}$  [$M_\odot$]')
    ax.legend(fontsize=legend_fs); ax.set_xlim(**xlim_t); ax.set_title('Cumulative ejecta mass')

    ax = axes[2]
    # Thin: instantaneous flux-weighted <Ye> crossing the sphere. Thick: <Ye>
    # of everything ejected so far -- the number to compare between models,
    # and only at equal t - t_merger, since it keeps rising as the wind arrives.
    for d, col, ls, name in _ejecta_lines(models):
        Yt, R = d['ye_avg'], d['rate']
        n = min(len(Yt), len(R))
        tt, yy, ww = Yt[:n, 0], Yt[:n, 1], R[:n, 1]
        dm = np.concatenate([[0.0], 0.5 * (ww[1:] + ww[:-1]) * np.diff(tt)])
        dmy = np.concatenate([[0.0], 0.5 * (ww[1:] * yy[1:] + ww[:-1] * yy[:-1]) * np.diff(tt)])
        M = np.cumsum(dm)
        ycum = np.where(M > 0, np.cumsum(dmy) / np.where(M > 0, M, 1.0), np.nan)
        ax.plot(tt, yy, color=col, ls=ls, lw=0.7, alpha=0.6)
        ax.plot(tt, ycum, color=col, ls=ls, lw=1.8,
                label=f'{name} cumulative' + ('' if multi else f' ({ycum[-1]:.3f})'))
    ax.plot([], [], color='gray', lw=0.7, alpha=0.6, label='instantaneous')
    ax.axhline(0.5, color='gray', ls=':', lw=0.8)
    ax.set_xlabel(xlab); ax.set_ylabel(r'$\langle Y_e \rangle_\mathrm{ej}$')
    ax.legend(fontsize=legend_fs); ax.set_xlim(**xlim_t); ax.set_ylim(0, 0.6); ax.set_title(r'Mass-weighted $\langle Y_e \rangle$ of the ejecta')

    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved {out}')


def plot_ye_theta(adir, radius):
    """2-D angular distribution of ejecta composition (Y_e vs theta).

    Reads the Mej_Ye_theta_{geo,Ber}.npy arrays; skips quietly if they are
    absent (analysis produced before this diagnostic was added).
    """
    try:
        theta = np.degrees(np.loadtxt(os.path.join(adir, 'theta_centers_2d.txt')))
        ye    = np.loadtxt(os.path.join(adir, 'Ye_centers.txt'))
        H_geo = np.load(os.path.join(adir, 'Mej_Ye_theta_geo.npy'))
        H_Ber = np.load(os.path.join(adir, 'Mej_Ye_theta_Ber.npy'))
    except (OSError, FileNotFoundError):
        print('  [ye_theta] Mej_Ye_theta_*.npy not found; skipping '
              '(re-run the ejecta analysis to generate it).')
        return

    fig, axes = plot_style.subplots(1, 2, panel_w=5.5, panel_h=4.2)
    fig.suptitle(f'BNS SPH — Ejecta composition vs angle  (r = {radius:g} M$_\\odot$)',
                 fontsize=13)
    vmax = max(H_geo.max(), H_Ber.max())
    norm = (LogNorm(vmin=vmax * 10**-LOG_DECADES, vmax=vmax) if vmax > 0 else None)
    for ax, H, title in ((axes[0], H_geo, 'geodesic'), (axes[1], H_Ber, 'Bernoulli')):
        pcm = ax.pcolormesh(theta, ye, np.ma.masked_less_equal(H, 0), cmap='viridis',
                            norm=norm, shading='auto')
        ax.set_xlabel(r'$\theta$ [deg]'); ax.set_ylabel(r'$Y_e$')
        ax.set_xlim(0, 180); ax.xaxis.set_major_locator(MultipleLocator(30))
        ax.set_ylim(0, 0.6); ax.set_title(f'{title}  ({_fmt_mass(H.sum())} $M_\\odot$)')
        fig.colorbar(pcm, ax=ax, label=r'$M_\mathrm{ej}$ per bin [$M_\odot$]')

    fig.tight_layout()
    out = os.path.join(adir, 'fig_ejecta_ye_theta.pdf')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved {out}')


def plot_neutrino(adir, t_ms, t_merger_ms, nu_xlabel, nu_xlim, radius):
    def load(f): return np.loadtxt(os.path.join(adir, f))

    SPECIES  = ['nue', 'nua', 'nux', 'anux']
    SP_LABEL = [r'$\nu_e$', r'$\bar\nu_e$', r'$\nu_x$', r'$\bar\nu_x$']
    SP_COLOR = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    Lnu  = {sp: load(f'Lnu_E_{sp}.txt') for sp in SPECIES}
    Ltot = load('Lnu_E_total.txt')
    Eav  = {sp: load(f'Eav_{sp}.txt')   for sp in SPECIES}

    # Two panels: (1) luminosity, log scale 1e51–3e54;  (2) mean energy, 0–60.
    fig, axes = plot_style.subplots(1, 2, panel_w=5.5, panel_h=4.2)
    fig.suptitle(f'BNS SPH — Neutrino emission  (r = {radius:g} M$_\\odot$)', fontsize=13)

    ax = axes[0]
    ax.semilogy(t_ms(Ltot[:, 0]), Ltot[:, 1], color='black', lw=1.5, label='total')
    for sp, lab, col in zip(SPECIES, SP_LABEL, SP_COLOR):
        d = Lnu[sp]; ax.semilogy(t_ms(d[:, 0]), d[:, 1], color=col, label=lab, lw=1.0)
    ax.axvline(t_merger_ms, color='gray', ls='--', lw=0.8)
    ax.set_xlabel(nu_xlabel); ax.set_ylabel(r'$L_{\nu,E}$  [erg s$^{-1}$]')
    ax.set_ylim(1e51, 3e54)
    ax.legend(fontsize=9, ncol=2); ax.set_xlim(**nu_xlim)
    ax.set_title('Energy luminosity')

    ax = axes[1]
    for sp, lab, col in zip(SPECIES, SP_LABEL, SP_COLOR):
        d = Eav[sp]; ax.plot(t_ms(d[:, 0]), d[:, 1], color=col, label=lab, lw=1.2)
    ax.axvline(t_merger_ms, color='gray', ls='--', lw=0.8)
    ax.set_xlabel(nu_xlabel); ax.set_ylabel(r'$\langle\varepsilon_\nu\rangle$  [MeV]')
    ax.set_ylim(0, 60)
    ax.legend(fontsize=9, ncol=2); ax.set_xlim(**nu_xlim)
    ax.set_title('Mean neutrino energy')

    fig.tight_layout()
    out = os.path.join(adir, 'fig_neutrino.pdf')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved {out}')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", required=True,
                        help="Directory holding the analysis .txt outputs.")
    parser.add_argument("--t-merger", type=float, required=True,
                        help="Merger time [M_sun] (from merger_time.txt).")
    parser.add_argument("--radius", type=float, default=300.0,
                        help="Neutrino SPH extraction radius (for plot titles). Default: 300.")
    parser.add_argument("--ej-radius", type=float, default=200.0,
                        help="Ejecta SPH extraction radius (for plot titles). Default: 200.")
    parser.add_argument("--from-merger", action="store_true",
                        help="Use t - t_merger [ms] as the neutrino x-axis (xlim left=0) "
                             "instead of absolute time.")
    parser.add_argument("--no-ejecta", action="store_true", help="Skip the ejecta figure.")
    parser.add_argument("--no-neutrino", action="store_true", help="Skip the neutrino figure.")
    args = parser.parse_args()

    adir = args.output_dir
    t_merger = args.t_merger

    if args.from_merger:
        def t_ms(t): return (t - t_merger) * MSUN_TO_MS
        nu_xlabel = r'$t - t_\mathrm{merger}$ [ms]'
        nu_xlim = {'left': 0}
        t_merger_ms = 0.0
    else:
        def t_ms(t): return t * MSUN_TO_MS
        nu_xlabel = r'$t$ [ms]'
        nu_xlim = {}
        t_merger_ms = t_merger * MSUN_TO_MS

    if not args.no_ejecta:
        model = [(None, None, load_ejecta_summary(adir, t_merger))]
        plot_ejecta(
            model, os.path.join(adir, 'fig_ejecta.pdf'),
            f'BNS SPH — Ejecta mass distributions  (r = {args.ej_radius:g} M$_\\odot$)')
        plot_ejecta_time(
            model, os.path.join(adir, 'fig_ejecta_time.pdf'),
            f'BNS SPH — Ejecta properties  (r = {args.ej_radius:g} M$_\\odot$, '
            f'$t_\\mathrm{{merger}}$ = {t_merger:g} M$_\\odot$)')
        plot_ye_theta(adir, args.ej_radius)
    if not args.no_neutrino:
        plot_neutrino(adir, t_ms, t_merger_ms, nu_xlabel, nu_xlim, args.radius)


if __name__ == "__main__":
    main()
