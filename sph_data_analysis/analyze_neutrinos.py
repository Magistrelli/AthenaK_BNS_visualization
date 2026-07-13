"""
Neutrino luminosity analysis from grouped rad_m1_moments or legacy
rad_m1_F / rad_m1_E / rad_m1_N spherical-surface VTK outputs.

Per snapshot, loads one grouped M1 file:
  bns.r=R.rad_m1_moments.NNNNN.vtk — E:s, N:s, Fx:s, Fy:s, Fz:s (all species)
  (older runs: bns.r=R.m1_mom.NNNNN.vtk, or split rad_m1_E/F/N files)
Legacy segments with three separate rad_m1_F / E / N files remain supported.
  bns.r=R.adm.NNNNN.vtk       — adm_gxx..gyz       (ADM 3-metric)
  bns.r=R.z4c_gauge.NNNNN.vtk — z4c_alpha, z4c_beta* (lapse and shift)
  (older runs: split z4c_alpha / z4c_betax / z4c_betay / z4c_betaz files)

R is either a single radius ('r=300.00') or, in newer output, a range
('r=100.00-500.00') holding many shells in one file; sph_vtk auto-detects
which style is on disk and extracts the shell nearest RADIUS.

Physics
-------
Energy luminosity per species s:
    L_E,s(t) = ∫ (α F^r_s - β^r E_s) · r² dΩ
             [code units → erg/s]
where F^i_s = γ^ij F_j,s, using the stored M1 flux covector components.
The radial projection uses the metric-normalized coordinate radial covector
used by THC: rhat_i = r_i / sqrt(g^ij_4 r_i r_j), where g^ij_4 is the
spatial block of the inverse 4-metric.
E_s and F_i,s are already densitized M1 evolved variables, so no extra √γ
factor is applied in the surface integral.
This matches the THC/AthenaK M1 calc_E_flux expression.
The angular area integral uses AthenaK's native SPH `weights` field:
    weights = r² dΩ  (Gauss-Legendre in μ = cos θ, uniform in φ)

Flux-weighted mean energy per species:
    <ε_s>(t) = ∫ F^r_E,s · r² dΩ / ∫ F^r_N,s · r² dΩ  → MeV
THC constructs the denominator from the M1 number flux on the detector sphere.
The AthenaK SPH files used here store E_s, F_i,s, and N_s but not the reconstructed
number-current direction fnu^a, so the post-process uses the same radial transport
velocity as the energy luminosity:
    F^r_N,s ≈ F^r_E,s · N_s / E_s.
This is exact when the energy and number currents are advected with the same
radial velocity, and is the appropriate free-streaming detector-sphere limit.

Species ordering (M1_TOTAL_NUM_SPECIES = 4):
    s = 0  νe   (electron neutrino)
    s = 1  ν̄e   (electron antineutrino)
    s = 2  νx   (heavy-flavour neutrino)
    s = 3  ν̄x   (heavy-flavour antineutrino)
Note: in bns_nurates, species 2 and 3 are symmetric (nux = anux).  The code
stores each separately; both are output here.

Required parfile output block (add for each desired radius):
    <outputXX>
    file_type  = sph
    variable   = rad_m1_moments
    dt         = 50.0
    radius     = 400.0
    ntheta     = 256
    id         = m1_mom

The three legacy rad_m1_F / rad_m1_E / rad_m1_N blocks are also supported.

Output files (written to OUTPUT_DIR):
    Lnu_E_{sp}.txt        time, energy luminosity per species [M_sun], [erg/s]
    Lnu_E_total.txt       time, sum over all 4 species         [M_sun], [erg/s]
    Eav_{sp}.txt          time, flux-weighted mean energy       [M_sun], [MeV]
    dEnu_dtheta_{sp}.txt  θ, time-integrated dE_ν/dθ           [rad],   [erg/rad]
"""

import os
import sys
import numpy as np
import inspect
from multiprocessing import cpu_count, get_all_start_methods, get_context

from snapshot_cache import SnapshotCache, TimeIndexCache, stamp as _src_stamp

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))
# plot-tools submodule lives at the repository root (one level above this folder).
_PLOTTOOLS = os.environ.get(
    "ATHPLOT_PATH", os.path.abspath(os.path.join(SCRIPTDIR, os.pardir, "plot-tools")))
sys.path.append(_PLOTTOOLS)
sys.path.append(SCRIPTDIR)  # utils_ejecta.py lives alongside this script

from sph_vtk import (SphericalData, set_default_radius, resolve_styles,
                     describe_styles, read_time, open_set,
                     M1_GROUPED_VARS)
from athplot.utils.units import conv_luminosity, conv_energy, cactus, cgs
from utils_ejecta import sum_over_time

# Workers must inherit the settings the bash wrapper assigns before main()
# (e.g. SPH_DIR_LIST), so use "fork" where the platform has it (newer
# Pythons no longer default to it); elsewhere, the default start method.
_MP = get_context("fork" if "fork" in get_all_start_methods() else None)


# ===================================================================
# Configuration  — edit this block before running
# ===================================================================

# --- Simulation ---
SIM_NAME = 'M135_M135_D45_DD2_VVLR_AMR_rererun'
SPH_DIR_BASENAMES = [
    'output_scratch/output-0000/sph',
    'output_scratch/output-0001/sph',
    'output_scratch/output-0002/sph',
]

# Extraction radius [code units, M_sun] — must match a radius in the parfile
RADIUS = 300.0

# --- Data directories ---
SPH_DIR_LIST = [os.path.join('../sim', SIM_NAME, b) for b in SPH_DIR_BASENAMES]

# --- Output ---
OUTPUT_DIR = os.path.join('output', SIM_NAME)

# --- Parallelism ---
# Number of worker processes for the snapshot loop.
# Set to 1 to disable multiprocessing (easier debugging).
N_WORKERS = min(8, cpu_count())

# Reuse per-snapshot results cached under OUTPUT_DIR/cache/ from earlier runs,
# and only process snapshots never seen before (see snapshot_cache.py). The
# reduction below always re-runs over cached + new results together, so outputs
# are identical either way. Set False (or SPH_CACHE=0 in run_analysis.sh) to
# force every snapshot to be recomputed.
USE_CACHE = True

# Process every STRIDE-th snapshot (1 = all).  Long runs can have thousands of
# snapshots; raise this so the loop stays tractable (luminosity curves stay
# smooth at STRIDE 5-10).  Set from run_analysis.sh.
STRIDE = 1


# ===================================================================
# Constants
# ===================================================================

# Species labels matching AthenaK's M1 ordering
SPECIES  = ['nue', 'nua', 'nux', 'anux']
NSPECIES = len(SPECIES)

# AthenaK stores N:s (neutrino number density) in EOS units (fm^-3) while
# E:s is in code energy-density units (M_sun^-2).  To form a consistent ratio,
# convert N from fm^-3 to code units (M_sun^-3):
#   1 code_length = G M_sun / c^2 = 1.477e18 fm
#   1 M_sun^3     = (1.477e18)^3 fm^3 = 3.222e54 fm^3
#   ∴ 1 fm^-3     = 3.222e54 M_sun^-3
_G_SI        = 6.67430e-11
_M_SUN_KG    = 1.98892e30
_C_SI        = 2.99792458e8
_FM_PER_M    = 1e15
N_CODE_PER_FM3 = (_G_SI * _M_SUN_KG / _C_SI**2 * _FM_PER_M)**3  # ≈ 3.222e54


# ===================================================================
# Main analysis
# ===================================================================

def _build_time_list(sph_dir, style, prefix):
    """Return sorted list of (time_float, path) for all matching VTK files."""
    result = []
    pat = os.path.basename(style.prefix(sph_dir)) + f'{prefix}.'
    for fname in os.listdir(sph_dir):
        if fname.startswith(pat) and fname.endswith('.vtk'):
            path = os.path.join(sph_dir, fname)
            # Read only the metadata line (avoids the full VTK parse)
            with open(path, 'r', encoding='latin_1') as f:
                f.readline()
                meta = f.readline()
            for tok in meta.split():
                if tok.startswith('time='):
                    result.append((float(tok[5:]), path))
                    break
    result.sort(key=lambda x: x[0])
    return result


def _read_vtk_time(path):
    """Header time of one spherical-surface VTK, or None. Reads 2 lines only."""
    try:
        with open(path, 'r', encoding='latin_1') as f:
            f.readline()
            meta = f.readline()
    except OSError:
        return None
    for tok in meta.split():
        if tok.startswith('time='):
            return float(tok[5:])
    return None


def _build_time_lists(styles, prefixes, cache=None):
    """Sorted (time, path) lists for several prefixes at once.

    One os.scandir per directory covers ALL prefixes -- the old code called
    _build_time_list once per prefix, i.e. five full listings of directories
    holding tens of thousands of files. Header reads are served from `cache`
    when the file is unchanged, which is what makes a resumed job cheap.
    """
    out = {p: [] for p in prefixes}
    for d, style in styles.items():
        # Both the radius tag and the variable grouping can differ between
        # segments, so rebuild the filename patterns for each directory.
        # Several logical fields can map to one file (z4c_gauge holds the
        # lapse and all three shift components), so one filename may feed
        # more than one entry of `out`.
        base = os.path.basename(style.prefix(d))
        pats = {}
        for p in prefixes:
            var = style.var_for(p)
            if var is not None:
                pats[p] = base + f'{var}.'
        try:
            entries = list(os.scandir(d))
        except OSError:
            continue
        for e in entries:
            name = e.name
            if not name.endswith('.vtk'):
                continue
            matched = [p for p in prefixes
                       if p in pats and name.startswith(pats[p])]
            if matched:
                t = None
                if cache is not None:
                    try:
                        st = e.stat()
                        t = cache.lookup(e.path, st.st_size, st.st_mtime_ns)
                    except OSError:
                        st = None
                    if t is None:
                        t = _read_vtk_time(e.path)
                        if t is not None and st is not None:
                            cache.put(e.path, st.st_size, st.st_mtime_ns, t)
                else:
                    t = _read_vtk_time(e.path)
                if t is not None:
                    for p in matched:
                        out[p].append((t, e.path))
    for p in prefixes:
        out[p].sort(key=lambda x: x[0])
    return out


def _nearest_path(time_list, t, tol=4.0):
    """Return path of snapshot with time nearest to t, or None if outside tol."""
    if not time_list:
        return None
    times = np.array([x[0] for x in time_list])
    idx   = np.searchsorted(times, t)
    # check both neighbours
    best_i, best_dt = None, np.inf
    for i in (idx - 1, idx):
        if 0 <= i < len(times):
            dt = abs(times[i] - t)
            if dt < best_dt:
                best_dt, best_i = dt, i
    if best_dt <= tol:
        return time_list[best_i][1]
    return None


def _inverse_spatial_metric(gxx, gyy, gzz, gxy, gxz, gyz):
    gamma_det = (gxx * gyy * gzz +
                 2 * gxy * gxz * gyz -
                 gxx * gyz**2 -
                 gyy * gxz**2 -
                 gzz * gxy**2)

    guxx = (gyy * gzz - gyz**2) / gamma_det
    guyy = (gxx * gzz - gxz**2) / gamma_det
    guzz = (gxx * gyy - gxy**2) / gamma_det
    guxy = (gxz * gyz - gxy * gzz) / gamma_det
    guxz = (gxy * gyz - gxz * gyy) / gamma_det
    guyz = (gxy * gxz - gxx * gyz) / gamma_det
    return guxx, guyy, guzz, guxy, guxz, guyz


def _raise_spatial_covector(Fx, Fy, Fz, gxx, gyy, gzz, gxy, gxz, gyz):
    guxx, guyy, guzz, guxy, guxz, guyz = _inverse_spatial_metric(
        gxx, gyy, gzz, gxy, gxz, gyz)

    Fux = guxx * Fx + guxy * Fy + guxz * Fz
    Fuy = guxy * Fx + guyy * Fy + guyz * Fz
    Fuz = guxz * Fx + guyz * Fy + guzz * Fz
    return Fux, Fuy, Fuz


def _metric_normalized_radial_projection(Vx, Vy, Vz,
                                         gxx, gyy, gzz, gxy, gxz, gyz,
                                         alpha, betax, betay, betaz,
                                         theta, phi, radius):
    """Project V^i along r_i / sqrt(g^ij_4 r_i r_j), matching THC."""
    guxx, guyy, guzz, guxy, guxz, guyz = _inverse_spatial_metric(
        gxx, gyy, gzz, gxy, gxz, gyz)
    oo_alpha2 = 1.0 / alpha**2

    # THC normalizes r_i with the inverse 4-metric.  Its spatial block is
    # g^ij_4 = gamma^ij - beta^i beta^j / alpha^2.
    g4uxx = guxx - betax * betax * oo_alpha2
    g4uyy = guyy - betay * betay * oo_alpha2
    g4uzz = guzz - betaz * betaz * oo_alpha2
    g4uxy = guxy - betax * betay * oo_alpha2
    g4uxz = guxz - betax * betaz * oo_alpha2
    g4uyz = guyz - betay * betaz * oo_alpha2

    rx = radius * np.sin(theta) * np.cos(phi)
    ry = radius * np.sin(theta) * np.sin(phi)
    rz = radius * np.cos(theta)

    rr_metric = np.sqrt(g4uxx * rx**2 + g4uyy * ry**2 + g4uzz * rz**2 +
                        2.0 * g4uxy * rx * ry +
                        2.0 * g4uxz * rx * rz +
                        2.0 * g4uyz * ry * rz)
    return (rx * Vx + ry * Vy + rz * Vz) / rr_metric


def _process_snapshot(args):
    """Per-snapshot worker — runs in a separate process.

    Imports are done inside the function so the worker is self-contained and
    can be pickled by multiprocessing without carrying module-level state.
    """
    (flux_path, ener_path, nden_path, adm_path, alpha_path,
     betax_path, betay_path, betaz_path,
     theta_edges) = args

    import numpy as _np
    from sph_vtk import open_set as _open_set
    from athplot.utils.units import conv_luminosity as _conv_lum, \
                                    conv_energy as _conv_e, \
                                    cactus as _cactus, cgs as _cgs

    # With the grouped output format E/N/F share one rad_m1_moments file and
    # the lapse and shift share one z4c_gauge file, so these paths are not all
    # distinct; open_set reads each file exactly once.
    _vtk = _open_set((flux_path, ener_path, nden_path, adm_path,
                      alpha_path, betax_path, betay_path, betaz_path))
    flux_vtk  = _vtk[flux_path]
    ener_vtk  = _vtk[ener_path]
    nden_vtk  = _vtk[nden_path]
    adm_vtk   = _vtk[adm_path]
    alpha_vtk = _vtk[alpha_path]
    betax_vtk = _vtk[betax_path]
    betay_vtk = _vtk[betay_path]
    betaz_vtk = _vtk[betaz_path]

    time  = flux_vtk.time
    radius = flux_vtk.r
    theta = flux_vtk.theta
    phi   = flux_vtk.phi
    surface_weights = flux_vtk.fields['weights']

    # -- Duplication correction --
    z4c_alpha_raw = alpha_vtk.fields['z4c_alpha']
    phi_zero    = (phi == phi.min())
    theta_north = (theta == theta.min())
    dup_corr = _np.ones(z4c_alpha_raw.shape)
    if z4c_alpha_raw[phi_zero & ~theta_north].max() > 1.0:
        dup_corr[phi_zero] = 0.5
    if z4c_alpha_raw[theta_north].max() > 1.0:
        dup_corr[theta_north] = 0.25

    gxx = adm_vtk.fields['adm_gxx'] * dup_corr
    gyy = adm_vtk.fields['adm_gyy'] * dup_corr
    gzz = adm_vtk.fields['adm_gzz'] * dup_corr
    gxy = adm_vtk.fields['adm_gxy'] * dup_corr
    gxz = adm_vtk.fields['adm_gxz'] * dup_corr
    gyz = adm_vtk.fields['adm_gyz'] * dup_corr
    z4c_alpha  = z4c_alpha_raw * dup_corr
    z4c_betax  = betax_vtk.fields['z4c_betax'] * dup_corr
    z4c_betay  = betay_vtk.fields['z4c_betay'] * dup_corr
    z4c_betaz  = betaz_vtk.fields['z4c_betaz'] * dup_corr

    nspecies = 4
    # N unit conversion: fm^-3 → M_sun^-3
    _G, _M, _C, _FM = 6.67430e-11, 1.98892e30, 2.99792458e8, 1e-15
    n_conv = (_G * _M / _C**2 / _FM)**3

    Lnu_E_list   = []
    Eav_list     = []
    dLdtheta_list = []

    dtheta = _np.diff(theta_edges)

    for s in range(nspecies):
        Fx = flux_vtk.fields[f'Fx:{s}'] * dup_corr
        Fy = flux_vtk.fields[f'Fy:{s}'] * dup_corr
        Fz = flux_vtk.fields[f'Fz:{s}'] * dup_corr
        E  = ener_vtk.fields[f'E:{s}']  * dup_corr
        N  = nden_vtk.fields[f'N:{s}']  * dup_corr

        Fux, Fuy, Fuz = _raise_spatial_covector(
            Fx, Fy, Fz, gxx, gyy, gzz, gxy, gxz, gyz)
        Fr            = _metric_normalized_radial_projection(
            Fux, Fuy, Fuz, gxx, gyy, gzz, gxy, gxz, gyz,
            z4c_alpha, z4c_betax, z4c_betay, z4c_betaz,
            theta, phi, radius)
        beta_r        = _metric_normalized_radial_projection(
            z4c_betax, z4c_betay, z4c_betaz,
            gxx, gyy, gzz, gxy, gxz, gyz,
            z4c_alpha, z4c_betax, z4c_betay, z4c_betaz,
            theta, phi, radius)
        radial_flux   = z4c_alpha * Fr - beta_r * E
        L_code        = _np.sum(radial_flux * surface_weights)
        L_erg         = _conv_lum(_cactus, _cgs, L_code)
        Lnu_E_list.append(L_erg)

        # THC's detector mean energy is energy luminosity divided by number
        # luminosity.  The exact M1 number flux uses fnu^a, which is not stored
        # in the current SPH output, so use the same radial transport velocity
        # as the energy flux: F_N^r ~= F_E^r * N/E.
        N_code = N * n_conv
        number_flux = _np.zeros_like(radial_flux)
        _np.divide(radial_flux * N_code, E, out=number_flux, where=(E > 0.0))
        Ndot_code = _np.sum(number_flux * surface_weights)
        eps_MeV = (_conv_e(_cactus, _cgs, L_code / Ndot_code) / (_cgs.eV * 1e6)
                   if Ndot_code > 0.0 else 0.0)
        Eav_list.append(eps_MeV)

        outgoing_L  = _np.where(radial_flux > 0.0, radial_flux, 0.0)
        dL_theta, _ = _np.histogram(
            theta.flatten(), bins=theta_edges,
            weights=outgoing_L.flatten() * surface_weights.flatten(),
        )
        dLdtheta_list.append(dL_theta / dtheta)

    return {'time': time, 'Lnu_E': Lnu_E_list,
            'Eav': Eav_list, 'dLdtheta': dLdtheta_list}


def _cached_worker_neutrino(pair):
    """(index, args) -> (index, result). Lets imap_unordered stay unordered
    while still telling us which snapshot each result belongs to."""
    idx, args = pair
    r = _process_snapshot(args)
    r['_src'] = _src_stamp(args[0])      # args[0] = this snapshot's flux VTK
    return idx, r


def main():
    # ------------------------------------------------------------------
    # Build snapshot index.  New segments store E/N/F together in m1_mom;
    # legacy segments use three rad_m1_* files.  Resolve the format per
    # directory so one analysis can span the transition.
    # ------------------------------------------------------------------
    # Auto-detect the sph output style: one radius per file
    # ('bns.r=300.00.m1_mom.NNNNN.vtk') or a whole range packed into one
    # ('bns.r=100.00-500.00.m1_mom.NNNNN.vtk'), from which the shell nearest
    # RADIUS is extracted.  Resolved per directory so a run that switched
    # format between segments still analyses end to end.
    set_default_radius(RADIUS)
    # var='rad_m1_F' only constrains which radius tag is eligible: it matches
    # any of rad_m1_F / m1_mom / rad_m1_moments, so a radius is considered
    # only if it actually has M1 output there.
    styles = resolve_styles(SPH_DIR_LIST, radius=RADIUS, var='rad_m1_F')
    print(describe_styles(styles, RADIUS))

    m1_index_map = {}
    for d, style in styles.items():
        names = os.listdir(d)
        # Scan legacy first so a grouped file wins if both exist for an index.
        for prefix in ("rad_m1_F",) + M1_GROUPED_VARS:
            pattern = os.path.basename(style.prefix(d)) + f"{prefix}."
            for fname in names:
                if fname.startswith(pattern) and fname.endswith(".vtk"):
                    idx = int(fname.split(".")[-2])
                    m1_index_map[idx] = (d, prefix)

    indices = sorted(m1_index_map)[:: max(1, STRIDE)]
    grouped_count = sum(m1_index_map[i][1] in M1_GROUPED_VARS
                        for i in indices)
    legacy_count = len(indices) - grouped_count
    print(f"Found {len(indices)} M1 snapshots at r={RADIUS} M_sun "
          f"({grouped_count} grouped, {legacy_count} legacy; STRIDE={STRIDE})")
    if not indices:
        raise RuntimeError(
            f"No rad_m1_moments / m1_mom / rad_m1_F VTK files found in "
            f"{SPH_DIR_LIST}.\n"
            "Add a rad_m1_moments sph output block (or the three legacy "
            "rad_m1_F / rad_m1_E / rad_m1_N blocks) to the parfile."
        )

    # Build sorted time lists for adm and z4c gauge outputs (counter may differ)
    print("Building time index for adm / z4c gauge files...")
    _tidx = TimeIndexCache(OUTPUT_DIR, f'timeindex-r{RADIUS:g}', enabled=USE_CACHE)
    _lists = _build_time_lists(
        styles,
        ('adm', 'z4c_alpha', 'z4c_betax', 'z4c_betay', 'z4c_betaz'), _tidx)
    _tidx.save()
    print(f"  [timeindex] {_tidx.hits} header read(s) reused, {_tidx.misses} read from disk")
    adm_times   = _lists['adm']
    alpha_times = _lists['z4c_alpha']
    betax_times = _lists['z4c_betax']
    betay_times = _lists['z4c_betay']
    betaz_times = _lists['z4c_betaz']
    print(f"  adm files:       {len(adm_times)}")
    print(f"  z4c_alpha files: {len(alpha_times)}")
    print(f"  z4c_beta  files: {len(betax_times)}")

    # ------------------------------------------------------------------
    # Build theta bins for dL/dtheta diagnostics.  Surface integrals use the
    # native AthenaK SPH weights stored in each VTK file.
    # ------------------------------------------------------------------
    first_dir, first_prefix = m1_index_map[indices[0]]
    first_flux_path = styles[first_dir].path(
        first_dir, first_prefix, indices[0])
    _fv = SphericalData(first_flux_path, radius=RADIUS)
    R_ACTUAL = float(_fv.r)
    if abs(R_ACTUAL - RADIUS) > 1e-6:
        print(f"  [radius] nearest available shell is r={R_ACTUAL:.4f} M_sun "
              f"(requested {RADIUS:g}); all integrals use r={R_ACTUAL:.4f}")
    theta_1d        = _fv.theta[0, :]
    theta_1d_sorted = np.sort(theta_1d)
    theta_edges     = np.concatenate([[0.0],
                                      0.5 * (theta_1d_sorted[:-1] + theta_1d_sorted[1:]),
                                      [np.pi]])
    theta_centers   = theta_1d_sorted
    del _fv

    # dup check on first snapshot with alpha
    _t0 = read_time(first_flux_path)
    _ap = _nearest_path(alpha_times, _t0)
    if _ap:
        _av = SphericalData(_ap)
        _a  = _av.fields['z4c_alpha']
        print(f"  [dup check] alpha range: min={_a.min():.4f}  max={_a.max():.4f}")

    # ------------------------------------------------------------------
    # Build worker argument list (pre-resolve adm/alpha paths by time)
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    # Incremental cache: reuse per-snapshot results computed by earlier runs.
    # This also skips the SERIAL SphericalData(flux_path).time header read
    # below, which on its own costs one file open per snapshot.
    # ------------------------------------------------------------------
    cache = SnapshotCache(
        OUTPUT_DIR, 'neutrino',
        params={'radius': RADIUS, 'r_actual': R_ACTUAL,
                'theta_edges': [float(x) for x in theta_edges],
                'species': list(SPECIES)},
        sources=(inspect.getsource(_process_snapshot),),
        enabled=USE_CACHE)
    cached = cache.load()

    results   = []
    jobs      = []
    skipped   = 0
    for index in indices:
        sph_dir, m1_prefix = m1_index_map[index]
        pfx = styles[sph_dir].prefix(sph_dir)
        flux_path = pfx + f"{m1_prefix}.{index:05d}.vtk"
        rec = cached.get(index)
        if rec is not None and rec.get('_src') == _src_stamp(flux_path):
            results.append(rec)          # note: skips the header read below
            continue
        t_flux = read_time(flux_path)   # 2-line header read, no data
        adm_path   = _nearest_path(adm_times,   t_flux)
        alpha_path = _nearest_path(alpha_times, t_flux)
        bx_path    = _nearest_path(betax_times, t_flux)
        by_path    = _nearest_path(betay_times, t_flux)
        bz_path    = _nearest_path(betaz_times, t_flux)
        if (adm_path is None or alpha_path is None or bx_path is None or
                by_path is None or bz_path is None):
            skipped += 1
            continue
        if m1_prefix in M1_GROUPED_VARS:
            ener_path = nden_path = flux_path
        else:
            ener_path = styles[sph_dir].path_for(sph_dir, 'rad_m1_E', index)
            nden_path = styles[sph_dir].path_for(sph_dir, 'rad_m1_N', index)
        jobs.append((index, (
            flux_path,
            ener_path,
            nden_path,
            adm_path, alpha_path,
            bx_path, by_path, bz_path,
            theta_edges,
        )))
    if skipped:
        print(f"  [skip] {skipped} snapshots with no adm/alpha/beta match")
    print(cache.describe(len(results), len(jobs)))

    # ------------------------------------------------------------------
    # Parallel snapshot loop (new snapshots only)
    # ------------------------------------------------------------------
    if jobs:
        print(f"Processing {len(jobs)} snapshots with N_WORKERS={N_WORKERS}...")
        if N_WORKERS > 1:
            with _MP.Pool(N_WORKERS) as pool:
                for idx, r in pool.imap_unordered(
                        _cached_worker_neutrino, jobs,
                        chunksize=max(1, len(jobs) // (N_WORKERS * 4))):
                    cache.add(idx, r)      # flushed every CHUNK -> a kill here
                    results.append(r)      # costs at most CHUNK snapshots
        else:
            for idx, args in jobs:
                _, r = _cached_worker_neutrino((idx, args))
                cache.add(idx, r)
                results.append(r)
        cache.flush()

    if not results:
        raise SystemExit('No snapshots to analyze.')

    # Progress summary
    results.sort(key=lambda x: x['time'])
    for i, res in enumerate(results):
        if i % 50 == 0:
            print(f"  snapshot {i:4d}/{len(results)}  t={res['time']:.2f} M_sun  "
                  f"L_nue={res['Lnu_E'][0]:.3e} erg/s")

    # ------------------------------------------------------------------
    # Unpack results into time-sorted accumulators
    # ------------------------------------------------------------------
    time_arr     = np.array([r['time']    for r in results])
    Lnu_E_snap   = [np.array([r['Lnu_E'][s]   for r in results]) for s in range(NSPECIES)]
    Eav_snap     = [np.array([r['Eav'][s]     for r in results]) for s in range(NSPECIES)]
    dLdtheta_snap = [np.array([r['dLdtheta'][s] for r in results]) for s in range(NSPECIES)]

    # Total luminosity across all species (already time-sorted)
    Lnu_total = sum(Lnu_E_snap[s] for s in range(NSPECIES))

    # Time-integrated dE_ν/dθ per species [erg/rad] via midpoint Riemann sum over time.
    # dLdtheta_snap is in code_lum/rad; time_arr is in M_sun (= code time);
    # their integral is in code_energy/rad → convert to erg/rad with conv_energy.
    _conv_e = conv_energy(cactus, cgs, 1.0)  # erg per code_energy ≈ 1.79e54
    dEnu_dtheta = [_conv_e * sum_over_time(dLdtheta_snap[s], time_arr)
                   for s in range(NSPECIES)]

    # ------------------------------------------------------------------
    # Print summary
    # ------------------------------------------------------------------
    print("\nPeak energy luminosities:")
    for s, sp in enumerate(SPECIES):
        print(f"  {sp:4s}: {Lnu_E_snap[s].max():.3e} erg/s")
    print(f"  total: {Lnu_total.max():.3e} erg/s")

    print("\nTime-averaged mean energies (last snapshot):")
    for s, sp in enumerate(SPECIES):
        print(f"  {sp:4s}: {Eav_snap[s][-1]:.2f} MeV")

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for s, sp in enumerate(SPECIES):
        # Energy luminosity time series
        np.savetxt(
            os.path.join(OUTPUT_DIR, f'Lnu_E_{sp}.txt'),
            np.column_stack((time_arr, Lnu_E_snap[s])),
            header=f'time[M_sun]    Lnu_E_{sp}[erg/s]',
            fmt='%.6e',
        )
        # Mean energy time series
        np.savetxt(
            os.path.join(OUTPUT_DIR, f'Eav_{sp}.txt'),
            np.column_stack((time_arr, Eav_snap[s])),
            header=(f'time[M_sun]    Eav_{sp}[MeV]  '
                    f'(flux-weighted mean energy at r={RADIUS})'),
            fmt='%.6e',
        )
        # Time-integrated angular distribution
        np.savetxt(
            os.path.join(OUTPUT_DIR, f'dEnu_dtheta_{sp}.txt'),
            np.column_stack((theta_centers, dEnu_dtheta[s])),
            header=f'theta[rad]    dEnu_dtheta_{sp}[erg/rad]  (time-integrated, outgoing only)',
            fmt='%.6e',
        )

    np.savetxt(
        os.path.join(OUTPUT_DIR, 'Lnu_E_total.txt'),
        np.column_stack((time_arr, Lnu_total)),
        header='time[M_sun]    Lnu_E_total[erg/s]  (sum over all 4 species)',
        fmt='%.6e',
    )

    print(f'\nAll outputs saved to {OUTPUT_DIR}')


if __name__ == '__main__':
    main()
