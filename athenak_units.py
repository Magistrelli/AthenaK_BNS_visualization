"""
Unit conversions for AthenaK BNS simulation outputs.

Modeled after scidata.units: explicit scale-factor constants, conv_* functions
for each physical quantity, and a plot_params() helper that returns
(scale, label, use_log, vmin, vmax) ready for use with BNS.py's render_frames.

AthenaK BNS uses GeometricKilometer code units (G=c=1, length unit = 1 km)
for geometry and radiation fields, with nuclear conventions for baryon/radiation
number densities (stored in fm^-3).  Hydro rest-mass density uses the standard
Cactus geometric convention (G=c=M_sun=1), which gives a slightly different
length unit (~1.477 km) but the same consistent density scale.

Supported unit systems
----------------------
'code'  : native AthenaK code units (default for all vars except density, which
          is currently plotted in CGS by BNS.py for historical reasons)
'cgs'   : Gaussian CGS (g, cm, s, erg)
'ngs'   : NGS = bns_nurates units (nm, g, s, MeV energy)

Conversion reference
--------------------
Opacity  κ       : code km^-1       →  NGS nm^-1       ÷ 1e12
                                    →  CGS cm^-1        ÷ 1e5
Rad. energy dens : code km^-2       →  NGS MeV/nm^3     × GK2NGS_ENEDENS
                                    →  CGS erg/cm^3      × CODE2CGS_ENEDENS
Rad. number dens : code fm^-3       →  NGS nm^-3        × 1e18
                                    →  CGS cm^-3        × 1e39
Number emiss. η₀ : code             →  NGS nm^-3 s^-1   × UNIT_ND_DOT
                                    →  CGS cm^-3 s^-1   × CODE2CGS_ND_DOT
Energy emiss. η₁ : code             →  NGS MeV nm^-3 s^-1 × UNIT_ED_DOT
                                    →  CGS erg cm^-3 s^-1  × CODE2CGS_ED_DOT
Rest-mass density: code (G=c=Msun=1)→  CGS g/cm^3       × CODE2CGS_DENS ≈ 6.178e17
Avg. ν energy E/N: code km^-2/fm^-3 →  MeV               × CODE2MEV_AVGENE
"""

import math

from matplotlib.colors import LogNorm, Normalize

# ---------------------------------------------------------------------------
# Physical constants (CGS, matching AthenaK unit_system.cpp)
# ---------------------------------------------------------------------------
c    = 2.99792458e10     # speed of light          [cm/s]
G    = 6.67408e-8        # gravitational constant   [cm^3 g^-1 s^-2]
MeV  = 1.6021766208e-6   # 1 MeV in erg
Msun = 1.98892e33        # solar mass               [g]

# ---------------------------------------------------------------------------
# GeometricKilometer (GK) unit scales  — G=c=1, length unit = 1 km
# ---------------------------------------------------------------------------
_GK_press  = G / c**4 * 1e10    # 1 erg/cm^3  expressed in km^-2
_GK_volume = 1e-15               # 1 cm^3      expressed in km^3
_GK_energy = G / c**4 * 1e-5    # 1 erg       expressed in km

# NGS (nm, g, s, MeV energy)
_NGS_volume = 1e21               # 1 cm^3  in nm^3   (1 cm = 1e7 nm)
_NGS_numden = 1e-21              # 1 cm^-3 in nm^-3

# Nuclear (fm, MeV)
_NUC_volume = 1e39               # 1 cm^3  in fm^3   (1 cm = 1e13 fm)

# ---------------------------------------------------------------------------
# Derived conversion constants
# ---------------------------------------------------------------------------

# Energy density: code (km^-2) → MeV/nm^3
GK2NGS_ENEDENS    = (_GK_volume / _NGS_volume) * (1.0/MeV / _GK_energy)

# Energy density: code (km^-2) → erg/cm^3
CODE2CGS_ENEDENS  = 1.0 / _GK_press

# Radiation number density: code (fm^-3) → nm^-3
CODE2NGS_NUMDENS  = 1e18         # 1 fm^-3 = 1e18 nm^-3

# Radiation number density: code (fm^-3) → cm^-3
CODE2CGS_NUMDENS  = 1e39         # 1 fm^-3 = 1e39 cm^-3

# Opacity: code (km^-1) → nm^-1
CODE2NGS_OPACITY  = 1e-12        # 1 km^-1 = 1e-12 nm^-1

# Opacity: code (km^-1) → cm^-1
CODE2CGS_OPACITY  = 1e-5         # 1 km^-1 = 1e-5  cm^-1

# Number emissivity η₀: code → nm^-3 s^-1
_c_km       = c * 1e-5           # c in km/s  ≈ 2.998e5 km/s
UNIT_ND_DOT = CODE2NGS_NUMDENS * _c_km   # ≈ 2.998e23

# Number emissivity η₀: code → cm^-3 s^-1
CODE2CGS_ND_DOT   = UNIT_ND_DOT * _NGS_numden   # ÷ 1e21

# Energy emissivity η₁: code → MeV nm^-3 s^-1
UNIT_ED_DOT       = GK2NGS_ENEDENS * _c_km       # ≈ 2.265e29

# Energy emissivity η₁: code → erg cm^-3 s^-1
CODE2CGS_ED_DOT   = UNIT_ED_DOT * MeV * _NGS_numden

# Rest-mass density: code (G=c=M_sun=1) → g/cm^3
# Uses the Cactus length unit G*Msun/c^2 ≈ 1.477 km, consistent with
# the scale_factor=6.178e17 already in BNS.py.
_unit_len_cgs     = G * Msun / c**2            # ≈ 1.477e6 cm
CODE2CGS_DENS     = Msun / _unit_len_cgs**3    # ≈ 6.178e17

# Baryon number density: code rest-mass density → fm^-3, n_b = rho / m_b.
# AthenaK's CompOSE EOS uses the table's neutron mass "mn" as m_b
# (eos_compose.cpp); 939.565 MeV for both DD2 and SFHo.
M_B_MEV           = 939.565346
_m_b_g            = M_B_MEV * MeV / c**2       # ≈ 1.675e-24 g
CODE2FM3_NB       = CODE2CGS_DENS / _m_b_g / _NUC_volume

# Neutrino degeneracy parameter y = n (n/J)^3 / (4 pi/(hc)^3), the quantity
# bns_nurates' CalculateThickParamsFromM1 maps to eta (distribution.hpp).
# Dimensionless; with n [fm^-3] and J [MeV fm^-3] (the m1ff output) divide by
# 4 pi/(hc)^3 in fm.  Constants as in bns_nurates constants.hpp.
HC_MEV_FM         = 4.135667696e-21 * 2.99792458e23   # h*c ≈ 1239.84 MeV fm
FOURPI_HC3_FM     = 4.0 * math.pi / HC_MEV_FM**3      # [MeV^-3 fm^-3]
NU_Y_FD_MAX       = 64.0 / 81.0   # y3: Fermi-Dirac limit (eta -> inf); fallback above

# Average neutrino energy E/N: code ratio (km^-2)/(fm^-3) → MeV
CODE2MEV_AVGENE   = GK2NGS_ENEDENS / CODE2NGS_NUMDENS  # ≈ 7.556e5

# Fluid-frame moments rad_J / rad_n are written in NUCLEAR units (MeV/fm^3 and
# fm^-3), not code units -- J is rescaled by ene_conv on output so that
# <e>=J/n is already in MeV (see derived_variables.cpp ~L1261-1269).  So the
# fluid energy density converts to ngs/cgs by the fm^3 volume factor only,
# NOT by GK2NGS_ENEDENS (which is for the code-unit lab E).
NUC2NGS_ENEDENS   = 1e18          # MeV/fm^3 → MeV/nm^3   (1 fm^3 = 1e-18 nm^3)
NUC2CGS_ENEDENS   = MeV * 1e39    # MeV/fm^3 → erg/cm^3   (≈ 1.602e33)

# Baryon number density: rest-mass density ρ (code) → n_b = ρ/m_b [fm^-3].
# m_b = 939.56535 MeV is the baryon (neutron) mass read from DD2.athtab
# (table_scalars "mn"); this matches the sim's own nb = rho/GetBaryonMass().
_MB_MEV     = 939.56535
_mb_gram    = _MB_MEV * MeV / c**2               # MeV/c^2 → g  (≈ 1.675e-24)
CODE2NB_FM3 = CODE2CGS_DENS / _mb_gram * 1e-39   # code ρ → baryon fm^-3  (≈ 369)

# Magnetic field: code → Gauss.
# AthenaK GRMHD uses the geometrized Gaussian (magnetic pressure = b^2/2)
# convention in the same G=c=M_sun=1 system as the rest-mass density, so the
# magnetic energy density b^2/2 (code) maps to erg/cm^3 via CODE2CGS_DENS*c^2.
# Hence B_Gauss = sqrt(8*pi*u_cgs) = b_code * sqrt(4*pi*CODE2CGS_DENS)*c.
CODE2CGS_BFIELD   = c * math.sqrt(4.0 * math.pi * CODE2CGS_DENS)  # ≈ 8.35e19 G

# Time: code (G=c=M_sun=1) → seconds / milliseconds
_unit_time_s      = _unit_len_cgs / c          # ≈ 4.927e-6 s
CODE2MS           = _unit_time_s * 1e3         # code time → ms

# ---------------------------------------------------------------------------
# Explicit conversion functions  (mirror scidata.units.conv_*)
# ---------------------------------------------------------------------------

def conv_dens(val, units='cgs'):
    """Rest-mass density (code G=c=Msun=1 units) → target units."""
    if units == 'code': return val
    if units == 'cgs':  return val * CODE2CGS_DENS
    raise ValueError(f"density: unsupported units {units!r} (use 'code' or 'cgs')")


def conv_enedens(val, units='cgs'):
    """Radiation energy density E (code km^-2) → target units."""
    if units == 'code': return val
    if units == 'ngs':  return val * GK2NGS_ENEDENS
    if units == 'cgs':  return val * CODE2CGS_ENEDENS
    raise ValueError(f"enedens: unsupported units {units!r}")


def conv_numdens(val, units='cgs'):
    """Radiation/baryon number density N (code fm^-3) → target units."""
    if units == 'code': return val
    if units == 'ngs':  return val * CODE2NGS_NUMDENS
    if units == 'cgs':  return val * CODE2CGS_NUMDENS
    raise ValueError(f"numdens: unsupported units {units!r}")


def conv_opacity(val, units='cgs'):
    """Absorption/scattering opacity κ (code km^-1) → target units."""
    if units == 'code': return val
    if units == 'ngs':  return val * CODE2NGS_OPACITY
    if units == 'cgs':  return val * CODE2CGS_OPACITY
    raise ValueError(f"opacity: unsupported units {units!r}")


def conv_emissivity_N(val, units='cgs'):
    """Number emissivity η₀ (code) → target units."""
    if units == 'code': return val
    if units == 'ngs':  return val * UNIT_ND_DOT
    if units == 'cgs':  return val * CODE2CGS_ND_DOT
    raise ValueError(f"emissivity_N: unsupported units {units!r}")


def conv_emissivity_E(val, units='cgs'):
    """Energy emissivity η₁ (code) → target units."""
    if units == 'code': return val
    if units == 'ngs':  return val * UNIT_ED_DOT
    if units == 'cgs':  return val * CODE2CGS_ED_DOT
    raise ValueError(f"emissivity_E: unsupported units {units!r}")


def conv_avg_energy(val, units='MeV'):
    """Average neutrino energy E/N (code) → target units."""
    if units == 'code': return val
    if units in ('ngs', 'MeV'): return val * CODE2MEV_AVGENE
    if units == 'cgs':          return val * CODE2MEV_AVGENE * MeV
    raise ValueError(f"avg_energy: unsupported units {units!r}")


def conv_bfield(val, units='cgs'):
    """Magnetic field magnitude |B| (code G=c=Msun=1 units) → target units."""
    if units == 'code':        return val
    if units in ('cgs', 'ngs'): return val * CODE2CGS_BFIELD   # Gauss
    raise ValueError(f"bfield: unsupported units {units!r}")


def conv_time(val, units='cgs'):
    """Coordinate time (code G=c=Msun=1 units) → target units."""
    if units == 'code':        return val
    if units in ('cgs', 'ngs'): return val * CODE2MS           # milliseconds
    raise ValueError(f"time: unsupported units {units!r}")


# ---------------------------------------------------------------------------
# Per-variable render parameters
# ---------------------------------------------------------------------------
# Original BNS.py hardcoded defaults (for reference / regression checks)
# These are the exact values that were in BNS.py before unit-system support
# was added.  Preserved here so future comparisons can verify nothing drifted.
#
#   var          scale_factor  norm                  label
#   dens         6.178e17      LogNorm(6e8, 2e15)    r"$\rho~[\mathrm{g\,cm^{-3}}]$"
#   temperature  1.0           vmin=0, vmax=70       r"$T~[\mathrm{MeV}]$"
#   s_00         1.0           vmin=0, vmax=0.5      r"$Y_e$"
#   E:0          1.0           LogNorm(1e-11, 1e-8)  "E:0"
#   N:0          1.0           LogNorm(1e-10, 1e-2)  "N:0"
#   abs_0:0      1.0           LogNorm(1e-18, 1e-9)  r"$\kappa_{\mathrm{abs},0}$"
#   abs_1:0      1.0           LogNorm(1e-18, 1e-9)  r"$\kappa_\mathrm{abs}$"
#   eta_0:0      1.0           LogNorm(1e-10, 1e-4)  r"$\eta_{0,0}$"
#   eta_1:0      1.0           LogNorm(1e-10, 1e-4)  r"$\eta_{1,0}$"
#   con_H        1.0           LogNorm(1e-9,  1e-4)  r"$|\mathcal{H}|$"
#   avg_energy   1.0           LogNorm(1e-4,  1e-3)  "Average energy (code units)"
#
# Code-unit (vmin, vmax) defaults for each variable.  The values for
# 'dens' are back-computed from the existing BNS.py defaults (scale 6.178e17,
# CGS range 6e8–2e15).  Radiation code ranges match the LogNorms above.
_CODE_VMIN = {
    'dens':      6e8    / CODE2CGS_DENS,   # ≈ 9.7e-10
    'E:0':       1e-11,
    'N:0':       1e-7,
    'N:2':       1e-7,
    'J:0':       1e-9,    # fluid-frame energy density (rad_J) [MeV/fm^3]
    'n:0':       1e-10,   # fluid-frame number density  (rad_n) [fm^-3]
    'abs_0:0':   1e-13,
    'abs_1:0':   1e-13,
    'eta_0:0':   1e-10,
    'eta_1:0':   1e-10,
    'con_H':     1e-9,
    'avg_energy':1e-4,
    'nu_energy': 30.0 / CODE2MEV_AVGENE,
    'bmag':      1e12   / CODE2CGS_BFIELD,  # ≈ 1.2e-8
    'Y_nu':      1e-5,
    # bns_nurates thick-distribution reconstruction (already in target units;
    # conversion is baked into the derived function, so scale = 1 for these).
    'nu_y':      0.0,     # dimensionless; FD ceiling 64/81=0.79, overshoot shown
    'nu_eta_t':  -20.0,   # degeneracy mu/T, UNCAPPED (sim clamps at 20)
    'nu_temp_t': 0.0,     # spectral temperature [MeV]
    'nu_avge':   0.0,     # fluid-frame mean energy J/n [MeV]
    'nu_Ynu':    0.0,     # neutrino fraction per baryon n_nu/n_b (linear, like Y_e)
    # heavy-lepton (nux) variants (species 2, 0.5 mu+tau weighting)
    'nux_y':     0.0,     # ceiling 0.5*64/81 = 0.395
    'nux_eta_t': -20.0,
    'nux_temp_t':0.0,
    'nux_avge':  0.0,
}
_CODE_VMAX = {
    'dens':      2e15   / CODE2CGS_DENS,   # ≈ 3.2e-3
    'E:0':       1e-8,
    'N:0':       1e-3,
    'N:2':       1e-3,
    'J:0':       1e0,     # [MeV/fm^3]
    'n:0':       1e-3,    # [fm^-3]
    'abs_0:0':   1e2,
    'abs_1:0':   1e2,
    'eta_0:0':   1e-4,
    'eta_1:0':   1e-4,
    'con_H':     1e-4,
    'avg_energy':1e-3,
    'nu_energy': 150.0 / CODE2MEV_AVGENE,
    'bmag':      1e15   / CODE2CGS_BFIELD,  # ≈ 1.2e-5
    'Y_nu':      1.0,
    'nu_y':      0.8,   # just above the fallback y3 = 64/81 ~ 0.790
    'nu_eta_t':  50.0,    # wider than the sim's 20 clamp to reveal true degeneracy
    'nu_temp_t': 60.0,    # [MeV], matches the matter-temperature colour scale
    'nu_avge':   50.0,    # [MeV]
    'nu_Ynu':    0.1,     # neutrino fraction per baryon (Y_e uses 0-0.5)
    'nux_y':     0.75,    # ~half the nue range (nux ceiling is halved)
    'nux_eta_t': 50.0,
    'nux_temp_t':60.0,
    'nux_avge':  50.0,
}

# Scale factors: (var, units) → multiplicative factor applied to code values
_SCALES = {
    # rest-mass density
    ('dens',     'cgs'): CODE2CGS_DENS,
    # radiation energy density
    ('E:0',      'ngs'): GK2NGS_ENEDENS,
    ('E:0',      'cgs'): CODE2CGS_ENEDENS,
    # radiation number density
    ('N:0',      'ngs'): CODE2NGS_NUMDENS,
    ('N:0',      'cgs'): CODE2CGS_NUMDENS,
    ('N:2',      'ngs'): CODE2NGS_NUMDENS,
    ('N:2',      'cgs'): CODE2CGS_NUMDENS,
    # fluid-frame radiation moments -- J in MeV/fm^3 (nuclear), n in fm^-3
    ('J:0',      'ngs'): NUC2NGS_ENEDENS,
    ('J:0',      'cgs'): NUC2CGS_ENEDENS,
    ('n:0',      'ngs'): CODE2NGS_NUMDENS,
    ('n:0',      'cgs'): CODE2CGS_NUMDENS,
    # absorption opacities
    ('abs_0:0',  'ngs'): CODE2NGS_OPACITY,
    ('abs_0:0',  'cgs'): CODE2CGS_OPACITY,
    ('abs_1:0',  'ngs'): CODE2NGS_OPACITY,
    ('abs_1:0',  'cgs'): CODE2CGS_OPACITY,
    # number emissivity
    ('eta_0:0',  'ngs'): UNIT_ND_DOT,
    ('eta_0:0',  'cgs'): CODE2CGS_ND_DOT,
    # energy emissivity
    ('eta_1:0',  'ngs'): UNIT_ED_DOT,
    ('eta_1:0',  'cgs'): CODE2CGS_ED_DOT,
    # average neutrino energy (ngs and cgs both yield MeV)
    ('avg_energy','ngs'): CODE2MEV_AVGENE,
    ('avg_energy','cgs'): CODE2MEV_AVGENE,
    ('nu_energy', 'ngs'): CODE2MEV_AVGENE,
    ('nu_energy', 'cgs'): CODE2MEV_AVGENE,
    # magnetic field magnitude (cgs and ngs both yield Gauss)
    ('bmag',     'cgs'): CODE2CGS_BFIELD,
    ('bmag',     'ngs'): CODE2CGS_BFIELD,
}

# Axis / colorbar labels: (var, units) → LaTeX string
_LABELS = {
    # density
    ('dens',     'code'): r"$\rho~[\mathrm{km}^{-2}]$",
    ('dens',     'cgs' ): r"$\rho~[\mathrm{g\,cm^{-3}}]$",
    # radiation energy density
    ('E:0',      'code'): r"$E_\nu^{(0)}~[\mathrm{km}^{-2}]$",
    ('E:0',      'ngs' ): r"$E_\nu^{(0)}~[\mathrm{MeV\,nm^{-3}}]$",
    ('E:0',      'cgs' ): r"$E_\nu^{(0)}~[\mathrm{erg\,cm^{-3}}]$",
    # radiation number density
    ('N:0',      'code'): r"$N_\nu^{(0)}~[\mathrm{fm}^{-3}]$",
    ('N:0',      'ngs' ): r"$N_\nu^{(0)}~[\mathrm{nm^{-3}}]$",
    ('N:0',      'cgs' ): r"$N_\nu^{(0)}~[\mathrm{cm^{-3}}]$",
    ('N:2',      'code'): r"$N_\nu^{(2)}~[\mathrm{fm}^{-3}]$",
    ('N:2',      'ngs' ): r"$N_\nu^{(2)}~[\mathrm{nm^{-3}}]$",
    ('N:2',      'cgs' ): r"$N_\nu^{(2)}~[\mathrm{cm^{-3}}]$",
    # fluid-frame radiation moments (comoving J, n); J native = MeV/fm^3
    ('J:0',      'code'): r"$J_\nu^{(0)}~[\mathrm{MeV\,fm^{-3}}]$",
    ('J:0',      'ngs' ): r"$J_\nu^{(0)}~[\mathrm{MeV\,nm^{-3}}]$",
    ('J:0',      'cgs' ): r"$J_\nu^{(0)}~[\mathrm{erg\,cm^{-3}}]$",
    ('n:0',      'code'): r"$n_\nu^{(0)}~[\mathrm{fm}^{-3}]$",
    ('n:0',      'ngs' ): r"$n_\nu^{(0)}~[\mathrm{nm^{-3}}]$",
    ('n:0',      'cgs' ): r"$n_\nu^{(0)}~[\mathrm{cm^{-3}}]$",
    # absorption opacities
    ('abs_0:0',  'code'): r"$\kappa_{a,0}^{(0)}~[\mathrm{km}^{-1}]$",
    ('abs_0:0',  'ngs' ): r"$\kappa_{a,0}^{(0)}~[\mathrm{nm}^{-1}]$",
    ('abs_0:0',  'cgs' ): r"$\kappa_{a,0}^{(0)}~[\mathrm{cm}^{-1}]$",
    ('abs_1:0',  'code'): r"$\kappa_{a,1}^{(0)}~[\mathrm{km}^{-1}]$",
    ('abs_1:0',  'ngs' ): r"$\kappa_{a,1}^{(0)}~[\mathrm{nm}^{-1}]$",
    ('abs_1:0',  'cgs' ): r"$\kappa_{a,1}^{(0)}~[\mathrm{cm}^{-1}]$",
    # number emissivity
    ('eta_0:0',  'code'): r"$\eta_0^{(0)}~[\mathrm{code}]$",
    ('eta_0:0',  'ngs' ): r"$\eta_0^{(0)}~[\mathrm{nm^{-3}\,s^{-1}}]$",
    ('eta_0:0',  'cgs' ): r"$\eta_0^{(0)}~[\mathrm{cm^{-3}\,s^{-1}}]$",
    # energy emissivity
    ('eta_1:0',  'code'): r"$\eta_1^{(0)}~[\mathrm{code}]$",
    ('eta_1:0',  'ngs' ): r"$\eta_1^{(0)}~[\mathrm{MeV\,nm^{-3}\,s^{-1}}]$",
    ('eta_1:0',  'cgs' ): r"$\eta_1^{(0)}~[\mathrm{erg\,cm^{-3}\,s^{-1}}]$",
    # average energy (report in MeV for both ngs and cgs)
    ('avg_energy','code'): r"$\langle\varepsilon_\nu\rangle~[\mathrm{km^{-2}/fm^{-3}}]$",
    ('avg_energy','ngs' ): r"$\langle\varepsilon_\nu\rangle~[\mathrm{MeV}]$",
    ('avg_energy','cgs' ): r"$\langle\varepsilon_\nu\rangle~[\mathrm{MeV}]$",
    ('nu_energy', 'code'): r"$\langle\varepsilon_\nu\rangle~[\mathrm{km^{-2}/fm^{-3}}]$",
    ('nu_energy', 'ngs' ): r"$\langle\varepsilon_\nu\rangle~[\mathrm{MeV}]$",
    ('nu_energy', 'cgs' ): r"$\langle\varepsilon_\nu\rangle~[\mathrm{MeV}]$",
    # dimensionless / already physical
    ('temperature','code'): r"$T~[\mathrm{MeV}]$",
    ('temperature','ngs' ): r"$T~[\mathrm{MeV}]$",
    ('temperature','cgs' ): r"$T~[\mathrm{MeV}]$",
    ('s_00',      'code'): r"$Y_e$",
    ('s_00',      'ngs' ): r"$Y_e$",
    ('s_00',      'cgs' ): r"$Y_e$",
    ('Y_nu',      'code'): r"$Y_\nu$",
    ('nu_y',      'code'): r"$y_\nu$",
    ('con_H',     'code'): r"$|\mathcal{H}|$",
    ('con_H',     'ngs' ): r"$|\mathcal{H}|$",
    ('con_H',     'cgs' ): r"$|\mathcal{H}|$",
    # magnetic field magnitude
    ('bmag',      'code'): r"$|B|~[\mathrm{code}]$",
    ('bmag',      'ngs' ): r"$|B|~[\mathrm{G}]$",
    ('bmag',      'cgs' ): r"$|B|~[\mathrm{G}]$",
    # bns_nurates thick-distribution reconstruction (n, J -> y, eta_t, temp_t)
    ('nu_y',      'code'): r"$y = n^4/J^3/(4\pi/(hc)^3)$",
    ('nu_y',      'ngs' ): r"$y = n^4/J^3/(4\pi/(hc)^3)$",
    ('nu_y',      'cgs' ): r"$y = n^4/J^3/(4\pi/(hc)^3)$",
    ('nu_eta_t',  'code'): r"$\eta_t = \mu/T$ (reconstructed)",
    ('nu_eta_t',  'ngs' ): r"$\eta_t = \mu/T$ (reconstructed)",
    ('nu_eta_t',  'cgs' ): r"$\eta_t = \mu/T$ (reconstructed)",
    ('nu_temp_t', 'code'): r"$T_\nu~[\mathrm{MeV}]$ (reconstructed)",
    ('nu_temp_t', 'ngs' ): r"$T_\nu~[\mathrm{MeV}]$ (reconstructed)",
    ('nu_temp_t', 'cgs' ): r"$T_\nu~[\mathrm{MeV}]$ (reconstructed)",
    # fluid-frame mean neutrino energy <e> = J/n
    ('nu_avge',   'code'): r"$\langle\varepsilon_\nu\rangle_\mathrm{fl}~[\mathrm{MeV}]$",
    ('nu_avge',   'ngs' ): r"$\langle\varepsilon_\nu\rangle_\mathrm{fl}~[\mathrm{MeV}]$",
    ('nu_avge',   'cgs' ): r"$\langle\varepsilon_\nu\rangle_\mathrm{fl}~[\mathrm{MeV}]$",
    # neutrino number fraction per baryon Y_nu = n_nu / n_b
    ('nu_Ynu',    'code'): r"$Y_\nu = n_\nu/n_b$",
    ('nu_Ynu',    'ngs' ): r"$Y_\nu = n_\nu/n_b$",
    ('nu_Ynu',    'cgs' ): r"$Y_\nu = n_\nu/n_b$",
    # heavy-lepton (nux) reconstruction variants
    ('nux_y',     'code'): r"$y_{\nu_x} = n^4/J^3/(4\pi/(hc)^3)$",
    ('nux_y',     'ngs' ): r"$y_{\nu_x} = n^4/J^3/(4\pi/(hc)^3)$",
    ('nux_y',     'cgs' ): r"$y_{\nu_x} = n^4/J^3/(4\pi/(hc)^3)$",
    ('nux_eta_t', 'code'): r"$\eta_{t,\nu_x} = \mu/T$ (reconstructed)",
    ('nux_eta_t', 'ngs' ): r"$\eta_{t,\nu_x} = \mu/T$ (reconstructed)",
    ('nux_eta_t', 'cgs' ): r"$\eta_{t,\nu_x} = \mu/T$ (reconstructed)",
    ('nux_temp_t','code'): r"$T_{\nu_x}~[\mathrm{MeV}]$ (reconstructed)",
    ('nux_temp_t','ngs' ): r"$T_{\nu_x}~[\mathrm{MeV}]$ (reconstructed)",
    ('nux_temp_t','cgs' ): r"$T_{\nu_x}~[\mathrm{MeV}]$ (reconstructed)",
    ('nux_avge',  'code'): r"$\langle\varepsilon_{\nu_x}\rangle_\mathrm{fl}~[\mathrm{MeV}]$",
    ('nux_avge',  'ngs' ): r"$\langle\varepsilon_{\nu_x}\rangle_\mathrm{fl}~[\mathrm{MeV}]$",
    ('nux_avge',  'cgs' ): r"$\langle\varepsilon_{\nu_x}\rangle_\mathrm{fl}~[\mathrm{MeV}]$",
}

# Variables that use a linear (not log) colorscale
_LINEAR_VARS = {'temperature', 's_00', 'nu_energy',
                'nu_y', 'nu_eta_t', 'nu_temp_t', 'nu_avge', 'nu_Ynu',
                'nux_y', 'nux_eta_t', 'nux_temp_t', 'nux_avge'}


def plot_params(var, units='code'):
    """Return render parameters for *var* in *units*.

    Parameters
    ----------
    var   : str — variable name as used in BNS.py (e.g. 'dens', 'abs_1:0')
    units : str — 'code', 'cgs', or 'ngs'

    Returns
    -------
    scale : float      — multiply code values by this before plotting
    label : str        — colorbar / axis label (LaTeX)
    use_log : bool     — True → LogNorm, False → linear Normalize
    vmin  : float      — suggested colorbar minimum (in target units)
    vmax  : float      — suggested colorbar maximum (in target units)
    """
    scale   = _SCALES.get((var, units), 1.0)
    label   = _LABELS.get((var, units), _LABELS.get((var, 'code'), var))
    use_log = var not in _LINEAR_VARS

    vmin_code = _CODE_VMIN.get(var)
    vmax_code = _CODE_VMAX.get(var)
    vmin = vmin_code * scale if vmin_code is not None else None
    vmax = vmax_code * scale if vmax_code is not None else None

    return scale, label, use_log, vmin, vmax


def make_norm(var, units='code', vmin=None, vmax=None):
    """Convenience: return a matplotlib Normalize or LogNorm for *var*.

    Supplied *vmin*/*vmax* override the defaults from plot_params().
    """
    _, _, use_log, default_vmin, default_vmax = plot_params(var, units)
    lo = vmin if vmin is not None else default_vmin
    hi = vmax if vmax is not None else default_vmax
    if use_log:
        return LogNorm(vmin=lo, vmax=hi)
    return Normalize(vmin=lo, vmax=hi)
