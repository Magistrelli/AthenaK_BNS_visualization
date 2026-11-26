# AthenaK BNS Visualization

Plotting and movie-generation tools for **binary neutron star (BNS)** merger
simulations run with [AthenaK](https://github.com/IAS-Astrophysics/athenak).

The main driver renders time-history plots and 2-D slice frames (density,
temperature, electron fraction, magnetic field magnitude |B|,
neutrino/radiation fields, …) from a simulation's `output-XXXX` restart
segments, then stitches the frames into movies with `ffmpeg`. It also writes
1-D time series such as the global max |B| and its location (`bmax`).

A companion driver, `BNS_3d.sh`, renders the 3-D snapshots (`bns.prim3D`) as
PyVista isosurfaces (optionally threaded with magnetic-field / velocity
streamlines) and Matplotlib orthoslices, including |B|, current (‖∇×B‖),
vorticity, and the Q-criterion.

## Contents

| File | Purpose |
|------|---------|
| `BNS_all.py` | Main visualization driver. Fully configurable via command-line flags (`--help`). File names are keyed on `--jobname` (default `bns`), so the same driver serves BNS/BHNS/BBH runs; each frame series is written to its own `Figs/<name>_xy/` subfolder, and sections whose input files are absent are skipped. |
| `BNS_all.sh` | Convenience wrapper: reads paths from `config.ini`, builds the argument list (incl. `JOBNAME`), runs `BNS_all.py`, then `make_movies.sh`. |
| `BNS_3d.sh` | 3-D driver: renders `bns.prim3D` snapshots as PyVista isosurfaces + Matplotlib orthoslices, then makes movies. Uses `plot_3d_data/`. |
| `BNS_grid.sh` | 2-D comparison driver: one figure of many slice panels — models × diagnostics, or one model over time. Uses `plot_2d_grid/`. |
| `plot_2d_grid/` | 2-D panel-grid package (diagnostic registry, frame selection, layout engine with shared colorbars). See its own `README.md`. |
| `plot_3d_data/` | 3-D rendering package (data loader, style/derived-field helpers, isosurface, volume & orthoslice scripts). See its own `README.md`. |
| `BNS_compare.py` / `BNS_compare.sh` | Comparison of runs: history (`.hst`) variables of two or more runs overlaid, and for two runs `[run A] [run B] [ratio A/B]` slice panels for opacity/radiation (and other) fields. |
| `make_movies.sh` | Builds one MP4/GIF per frame subfolder (`Figs/<name>_xy/`) into `Figs/all_movies/`. |
| `athenak_units.py` | Unit conversions (code / CGS / NGS) and per-quantity plot parameters. |
| `plot_style.py` | Central Matplotlib house style (single source). `plot_style.use()` / `plot_style.subplots(nrows, ncols)` apply the fixed look and scale the figure size to the subplot grid (size items constant per panel). |
| `plot_time_per_cycle.py` | Standalone tool: wall-clock time-per-cycle distribution from `.out` logs. |
| `sph_data_analysis/` | SPH ejecta + neutrino analysis (spherical-surface VTK output). See its section below. |
| `config.example.ini` | Template for machine-specific paths → copy to `config.ini`. |
| `plot-tools/` | Git submodule: the [`athplot`](https://github.com/jfields7/plot-tools) package used to read AthenaK binary/VTK output. |

## Setup

```bash
# 1. Clone with the submodule
git clone --recurse-submodules <this-repo-url>
cd AthenaK_BNS_visualization
# (or, if already cloned:  git submodule update --init --recursive)

# 2. Install Python dependencies
pip install -r requirements.txt
pip install -e plot-tools          # installs athplot + its deps (vtk, PyCompOSE)

# 3. Create your local config (paths + environment — the ONE machine file)
cp config.example.ini config.ini
$EDITOR config.ini                 # set env_setup and simpath (figpath optional)
```

`ffmpeg` is required for movie generation. All machine-specific settings —
including how the environment is activated — live in `config.ini`; you do **not**
edit the `.sh` scripts to switch machines. Set the `env_setup` key to whatever
loads your toolchain (`module load python`, `conda activate …`, `source …/activate`);
the scripts run it via the shared `load_config.sh` helper.

**3-D pipeline (`BNS_3d.sh`)** additionally needs `pyvista`, `h5py`, `vtk`, and
`imageio-ffmpeg` (`pip install -r plot_3d_data/requirements.txt imageio-ffmpeg`).
On Aurora the `python` module is a venv (no `--user` installs), so install these
into a directory and point `pythonpath_extra` at it in `config.ini`:

```bash
module load python/3.12.12 py-numpy py-matplotlib py-scipy py-pip
python3.12 -m pip install --target ~/.python-extra/py312 pyvista imageio-ffmpeg h5py
# then in config.ini:  pythonpath_extra = /home/<you>/.python-extra/py312
```

## Usage

### Full pipeline (frames + movies)

```bash
bash BNS_all.sh
```

This reads `simpath`/`figpath` from `config.ini`, renders the `SECTIONS` selected
inside `BNS_all.sh`, writes frames to `figpath` (default `<simpath>/Figs`), and
assembles movies into `<figpath>/all_movies/`.

### Several simulations at once

List each simulation as its own `simpath` line in `config.ini`, optionally with a
legend label after `|` (default: the folder name); commented lines are skipped:

```ini
simpath = /path/to/run_A | base
simpath = /path/to/run_B | nucl
#simpath = /path/to/run_C
```

`BNS_all.sh`, `BNS_3d.sh` and `sph_data_analysis/run_analysis.sh` then run once per
simulation, one after the other, and report the ones that failed; the per-run
keys `figpath` and `t_merger` must stay blank. `BNS_compare.sh` and
`sph_data_analysis/run_comparison.sh` compare these simulations unless their own
`RUNS`/`MODELS` lists or command-line arguments say otherwise.

```bash
bash run_all.sh                    # all steps below, for all simulations
bash run_all.sh --only analysis,run_comparison
```

`run_all.sh` runs, in order: `BNS_all.sh` and `run_analysis.sh` for each
simulation, then `BNS_compare.sh` and `run_comparison.sh` once for all of them
(the latter only for simulations with an analysis directory), and ends with a
summary of what ran, was skipped or failed.

### Driver directly

`BNS_all.py` takes all configuration on the command line — no editing required:

```bash
python3 BNS_all.py --simpath /path/to/run \
                   --figpath /path/to/run/Figs \
                   --units cgs \
                   --sections history rho_max nu_energy \
                   --no-trackers --skip-existing
```

Run `python3 BNS_all.py --help` for the full list of sections and flags.

### Output-layout compatibility

No format switch is needed: readers detect the files in each `output-XXXX`
segment, so a run may change layouts between restarts.

| Data | Older split outputs | Grouped outputs |
|------|---------------------|-----------------|
| Radiation slices | `rad_E_xy`, `rad_N_xy` (also legacy `rad_m1_E/N`) | `m1_xy` or `m1_mom_xy` |
| Opacity slices | `eta_0_xy`, `eta_1_xy`, `abs_0_xy`, `abs_1_xy` | `m1opac_xy` or `m1_opac_xy` |
| Fluid-frame slices | `rad_J_xy`, `rad_n_xy` | `m1ff_xy` |
| SPH radiation | `rad_m1_E`, `rad_m1_F`, `rad_m1_N` | `rad_m1_moments` or `m1_mom` |
| SPH gauge | `z4c_alpha`, `z4c_betax/y/z` | `z4c_gauge` |

Slice names also apply to the corresponding `xz`/`yz` planes supported by
each driver. `BNS_compare.py` and the panel-grid driver retain their existing
xz-only fallback for the oldest unqualified binary IDs.
`BNS_all.py`, comparisons, and panel grids all recognize `m1_xy`/`m1opac_xy`
as written by `LR_dU.par`. Grouped files are preferred when both slice layouts
are present. Fluid-frame derived plots use `J`/`n`, not lab-frame `E`/`N`.

SPH analysis accepts both `r=300.00` and ranges such as `r=10.00-500.00`,
selecting the nearest stored shell to the configured extraction radius.
Ejecta integration uses each snapshot's angular grid, including when its
resolution changes between segments. Grouped `rad_m1_moments` already contains
E, F, and N: separate output blocks are not required. Quantities absent from
both layouts still require enabling the corresponding simulation outputs.

Run the compatibility regression tests in the configured Python environment:

```bash
python3 -m unittest discover -s tests -v
```

### 2-D comparison grids

```bash
bash BNS_grid.sh
```

Builds **one** figure holding many 2-D slice panels, so several models and
diagnostics can be compared side by side instead of flipping between PNGs.
`ROWS` and `COLS` each take one of `model`, `diagnostic` or `time`, and whatever
is left over is swept into one figure per value — which covers models ×
diagnostics, one model over time, and one diagnostic across models, with a
single mechanism.

Sample with `TIMES=(4000)` to compare runs at the same physical time (each
model uses its nearest frame, and the resolved frame is printed per model), or
`FRAMES=(1000)` to pick by index — index N is exact but run-relative, so runs
with different output cadences land at different times.

Colour scales are fixed per diagnostic (from `athenak_units.plot_params`,
overridable via `CLIM[density]="1e10,1e15"`), so every panel of a diagnostic is
directly comparable and one shared colorbar is correct. `COLORBAR=auto` places
it where the geometry allows: under each column when `COLS=diagnostic`, right of
each row when `ROWS=diagnostic`, once for the figure otherwise. `SHAREX`/`SHAREY` (on by default) drop both the
axis label and the tick labels on interior panels, so the bottom row and left
column carry them — one x label per column, not one caption for the figure.
Panel geometry is explicit rather than `constrained_layout`, so `PANEL_SIZE` and
`SPACING` (gap as a fraction of panel size) are exact, and each colorbar is
built from its panels' positions so it spans them precisely. `GRID=false`
disables the faint axis grid lines. `DRY_RUN=true bash BNS_grid.sh` prints the
command without rendering.

### 3-D isosurfaces & orthoslices

```bash
bash BNS_3d.sh
```

Finds every `bns.prim3D.*.bin` under `<simpath>/output-*/`, renders the `PLOTS`
(`orthoslice`, `isosurface`, `isosurface_colorT`, `isosurface_colorYe`,
`volume`) for each `SCALARS` entry (`log_dens`, `temperature`,
`log_bmag`, and optionally `current`, `vorticity`, `qcriterion`), writes PNGs to
`<simpath>/Figs_3d/`, and stitches movies. Edit `SCALARS`, `STREAMLINES`
(`"B"`/`"velocity"` field lines on isosurfaces), `STRIDE`, `JOBS`, `UNITS`,
`RESOLUTION`/`CROP`, and the per-scalar `ISO_LEVELS`/`ISO_CLIM` at the top of the
script. Parallelism uses `xargs -P` (no GNU `parallel` needed) and movies use the
`imageio-ffmpeg` binary if no system `ffmpeg` is present.

Run `DRY_RUN=true bash BNS_3d.sh` to print the render commands without spending
node time on them.

For smooth hydrodynamic fields (`rho`, temperature, and `Y_e`), the driver uses
linear AMR reconstruction, 320 cells across the longest crop axis, SSAA, smooth
surface normals, and volume-preserving Taubin contour smoothing. Set
`INTERPOLATION=nearest` and `SMOOTH_ITERATIONS=0` when exact native cell edges
are more important than presentation. The 320^3 preset is substantially more
memory-intensive than 256^3, so the default `JOBS` is limited to 4.

**Composition and volume modes.** `isosurface_colorYe` paints Y_e on
`log10(rho) = 11, 13` density surfaces (`YE_LEVELS`/`YE_CLIM`/`YE_CMAP`), and
`volume` renders density as a volume (`VOL_*`). The split is deliberate: a
volume rendering cannot show composition once the ejecta is optically thick,
since every ray then exits through tens of M of `Y_e ~ 0.4` outflow and the
frame saturates. Density survives volume rendering only because its colour ramp
paints tenuous material dark. Y_e is read from `ye` or from the `s_00` passive
scalar, whichever the snapshot carries.

Volume frames are the memory-heavy ones (a 480³ float32 field is ~422 MiB before
VTK's own copy), so `BNS_3d.sh` caps their concurrency at `VOL_MAX_JOBS`
independently of `JOBS`.

**Merger-relative timing and annotations.** Set `T0` to the merger time in code
units and every frame is labelled `t_merger +X.XXX ms`. `FONT` picks the
annotation typeface (`sans`, `serif`, `cm`, or VTK's `arial`/`times`/`courier`),
and `ANNOTATE=mathtext` (the default) composites the clock, colour bar and ruler
with Matplotlib so the maths is typeset properly — VTK's own text cannot set
`log₁₀ρ` or `Y_e` and would print `log10 rho [g cm^-3]` instead. `SCALEBAR="auto"`
adds a self-relabelling 1/2/5 ruler, `AXES=true` an orientation triad, and
`CAMERA_TRACK` interpolates a zoom across the sequence. The ruler is exact
because it forces an orthographic camera, where `PARALLEL_SCALE` is literally
the half-height of the view in code units.

The `ISO_LEVELS`/`ISO_CLIM` defaults assume `UNITS=cgs` (|B| in Gauss, ρ in
g cm⁻³); switch the commented code-unit ranges if you set `UNITS=code`. The
underlying scripts can also be called directly — see `plot_3d_data/README.md`.

### Comparing runs

```bash
bash BNS_compare.sh          # or: python3 BNS_compare.py --run A --run B ...
```

Overlays the history (`.hst`) variables of all runs in one figure. With exactly
two runs it also renders `[run A] [run B] [ratio A/B]` slice panels for matched
output times, rasterizing each run's (possibly different) AMR grid onto a
common grid first.

### Movies only

```bash
bash make_movies.sh <fps> <stride> <FIGPATH>
# e.g. every 10th frame at 15 fps:
bash make_movies.sh 15 10 /path/to/run/Figs
```

### Reruns only redo what changed

With `SKIP_EXISTING=true` (default in `BNS_all.sh` and `BNS_compare.sh`) a frame
is only rendered again if its input files (name, size, modification time) or
the plot settings of its series changed; each frame folder records this in
`frames.json` (`frame_cache.py`). A movie is only re-encoded when its frames
or movie settings changed (`movie_tools.sh`, shared by `make_movies.sh` and
`BNS_compare.sh`). Afterwards the frames are packed, unchanged, into
`<folder>/frames.tar` (`PACK_FRAMES=1`; `tar -xf frames.tar` restores the PNGs),
and unpacked again automatically when new frames need a new movie. The SPH
analysis (`run_analysis.sh`) likewise reuses its per-snapshot cache.

### Time-per-cycle diagnostic

```bash
python3 plot_time_per_cycle.py --sim-dir /path/to/run \
                               --sim-name my_run \
                               --output-dirs output-0000 output-0001 \
                               --no-show
```

## Configuration

All **machine-specific** settings live in a single file, `config.ini` — copied
from `config.example.ini` and **gitignored**, so nothing machine-local gets
committed. It holds both:

* **paths** — `simpath` (one line per simulation, see above), optional
  `figpath`, and (for `sph_data_analysis`) `eos_table` / `pythonpath_extra`, and
* **environment** — `env_setup`, the one command run to activate your toolchain
  (`module load python`, `conda activate …`, `source …/activate`).

Every wrapper (`BNS_all.sh`, `BNS_3d.sh`, `BNS_compare.sh`,
`sph_data_analysis/run_analysis.sh`) sources the shared `load_config.sh` helper,
which parses `config.ini` (values are read as data; only `env_setup` is executed)
and activates the environment. So to move to a new machine you edit `config.ini`
only — never the scripts. They exit with an error if `config.ini` is missing.

All **run settings** — which `SECTIONS` to plot, `UNITS`, tracker/domain flags,
and movie `FPS`/`STRIDE` — stay in `BNS_all.sh`, where you edit them per run.

The plot-tools location is resolved automatically (the bundled submodule), or
can be overridden with the `ATHPLOT_PATH` environment variable.

## SPH data analysis (`sph_data_analysis/`)

Post-processing of the spherical-surface (`sph`) VTK output: dynamical **ejecta**
properties and **neutrino** emission, plus a **merger-time** estimate from the
gravitational waveform.

| File | Purpose |
|------|---------|
| `run_analysis.sh` | Orchestrator: reads paths from `../config.ini`, then runs the merger-time, ejecta, neutrino, plotting, and accretion steps. |
| `merger_time.py` | Merger time = peak of the (2,2) strain amplitude \|h₂₂\|, with h obtained from psi4₂₂ by fixed-frequency integration, merged across all `output-*/waveforms`. A restart that stops before merger falls back to its restart parent's waveform (`output-0000/parent`); `t_merger` in config.ini overrides. Writes `merger_time.txt`. |
| `analyze_ejecta.py` | Per-snapshot ejecta diagnostics: mass-ejection rate, and mass distributions vs v∞, temperature, Yₑ, entropy-per-baryon (from EOS Q2), θ, plus the 2-D Yₑ–θ composition map. Needs the EOS table. |
| `analyze_neutrinos.py` | Neutrino energy luminosity, mean energy, and angular distribution from grouped `rad_m1_moments` or split `rad_m1_E/F/N`. |
| `analyze_accretion.py` | Post-merger accretion rate from baryon bookkeeping (`Ṁ_accr = −dM_b/dt + Ṁ_ej(t−r/⟨v⟩)`), reusing the `.hst` files + ejecta outputs. Writes `fig_accretion.png` / `Mdot_accr_*.txt`. |
| `plot_analysis.py` | Builds `fig_ejecta.pdf` (1×5 mass distributions vs v∞/T/Yₑ/entropy, log Mₑⱼ 1e-6–1e-2 M⊙, and dM/dθ), `fig_ejecta_time.pdf` (1×3 histories: rate, cumulative mass, instantaneous/cumulative ⟨Yₑ⟩), `fig_ejecta_ye_theta.pdf` (Yₑ vs θ), and `fig_neutrino.pdf` (log luminosity 1e51–3e54, mean energy 0–60). |
| `run_comparison.sh` | Overlays the analysis outputs of 2–4 runs: neutrino luminosity and mean energy and the ejecta figures `fig_ejecta.pdf` and `fig_ejecta_time.pdf` (all via `comparison_plot.py`), with times relative to each run's merger. |
| `utils_ejecta.py` | Shared spherical-integration helpers. |

```bash
cd sph_data_analysis
bash run_analysis.sh              # merger time + ejecta + neutrino + plots + accretion
bash run_analysis.sh --neutrino --plot   # selected steps only
bash run_analysis.sh --accretion         # accretion only (needs prior ejecta outputs)
```

Defaults: a **single extraction radius `r = 300`** for both ejecta and neutrino
analysis, and **all** `output-*/sph` segments are processed (sorted, like
`BNS_all.sh`). Paths (`simpath`, `eos_table`, optional `pythonpath_extra`) come
from `config.ini`; run settings (radius, `N_WORKERS`, post-merger window, etc.)
are edited at the top of `run_analysis.sh`. Outputs are written to
`<simpath>/analysis/`.

Requires `h5py` (for the EOS table) and `vtk` (via `athplot.load_sph_vtk`) in
addition to the base dependencies. `vtk` is usually **not** in the bare
module/conda Python, so point `env_setup` at an environment that has it (or add
it via `pythonpath_extra`). On Aurora: `module load py-h5py` + `pythonpath_extra`.
On Perlmutter: a venv built from the `python` module with `pip install vtk`
(then `env_setup = source <venv>/bin/activate`) provides vtk + h5py + the rest.

## Notes

- **Plot styling** is defined once in `plot_style.py`; edit it there to change
  fonts/ticks/line widths everywhere. Fixed style items (serif font, inward
  minor ticks, frameless legend) and the per-panel size items are constant;
  only the figure size scales with the subplot grid, via
  `plot_style.figsize(nrows, ncols)` / `plot_style.subplots(nrows, ncols)`.
  Real LaTeX text is off by default (`plot_style.use(usetex=True)` to enable).
- `core.*` dumps, `__pycache__/`, and generated output (`Figs/`, `*.png`,
  `*.mp4`, …) are gitignored.
- This repo is focused on the visualization scripts. The separate
  [`LargeSim`](https://github.com/computationalrelativity/LargeSim) workflow
  repo and the older Perlmutter script variants are intentionally excluded.
