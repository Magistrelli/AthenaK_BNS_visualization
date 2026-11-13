"""Regression tests for split/grouped slices and single/multi-radius SPH."""

import contextlib
import io
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from multiprocessing.dummy import Pool
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "plot_2d_grid"),
                str(ROOT / "sph_data_analysis")]

import BNS_all as bns
import BNS_compare as compare
import diagnostics
import griddata
import analyze_ejecta as ejecta
import analyze_neutrinos as neutrinos
import sph_vtk as sph
import sph_extract


def write_bin(path, time, fields, plane):
    path.parent.mkdir(parents=True, exist_ok=True)
    par = (b"<mesh>\nnghost=0\nx1min=-1\nx1max=1\nx2min=-1\n"
           b"x2max=1\nx3min=-1\nx3max=1\n")
    indices = (0, 1, 0, 1, 0, 0) if plane == "xy" else (0, 1, 0, 0, 0, 1)
    with path.open("wb") as f:
        f.write(("Athena binary output version=1.1\n  size of preheader=5\n"
                 f"  time={time}\n  cycle=1\n  size of location=8\n"
                 f"  size of variable=4\n  number of variables={len(fields)}\n"
                 f"  variables: {' '.join(fields)}\n"
                 f"  header offset={len(par)}\n").encode())
        f.write(par)
        f.write(struct.pack("@6i", *indices))
        f.write(struct.pack("@4i", 0, 0, 0, 0))
        f.write(struct.pack("=6d", -1, 1, -1, 1, -1, 1))
        for value in fields.values():
            f.write(np.full(4, value, dtype="=f4").tobytes())


class SliceLayouts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.pool = patch.object(bns.multiprocessing, "Pool", Pool)
        self.pool.start()
        self.addCleanup(self.pool.stop)
        compare._FRAME_CACHE.clear()
        self.dirs = []
        for idx, style in enumerate(("legacy", "split", "grouped", "current")):
            d = self.root / f"output-{idx:04d}"
            self.dirs.append(str(d))
            for plane in ("xy", "xz"):
                fields = {"E:0": 2e-9, "N:0": 1e-5, "N:2": 1e-5,
                          "eta_0:0": 1e-8, "eta_1:0": 1e-8,
                          "abs_0:0": 1e-12, "abs_1:0": 1e-12,
                          "J:0": 2e-4, "n:0": 1e-5}
                groups = {}
                for var, value in fields.items():
                    base = var.split(":")[0]
                    if base in ("J", "n"):
                        prefix = "m1ff" if style == "current" else f"rad_{base}"
                    elif style in ("legacy", "split"):
                        prefix = f"rad_{base}" if base in ("E", "N") else base
                    elif base in ("E", "N"):
                        prefix = "m1" if style == "current" else "m1_mom"
                    else:
                        prefix = "m1opac" if style == "current" else "m1_opac"
                    prefix = f"{prefix}_{plane}"
                    if style == "legacy" and base not in ("J", "n"):
                        prefix = f"rad_m1_{base}"
                    groups.setdefault(prefix, {})[var] = value
                groups[f"prim_{plane}"] = {"dens": 1e-4}
                for prefix, variables in groups.items():
                    # The oldest unqualified binary IDs have one slice only.
                    if style == "legacy" and plane == "xz" and prefix.startswith("rad_m1_"):
                        continue
                    write_bin(d / "bin" / f"bns.{prefix}.00001.bin",
                              idx + 1, variables, plane)

    def test_main_comparison_and_grid_detect_each_segment(self):
        names = ("rad_E", "rad_N", "rad_N2", "rad_abs", "rad_abs0",
                 "rad_eta0", "rad_eta1")
        argv = ["BNS_all.py", "--simpath", str(self.root), "--sections", *names]
        with patch.object(sys, "argv", argv), \
                patch.object(bns, "render_frames_all") as render, \
                patch.object(bns, "load_trackers_all", return_value=[]), \
                contextlib.redirect_stdout(io.StringIO()):
            bns.main()
        self.assertEqual(render.call_count, len(names))
        for call in render.call_args_list:
            frames = bns._collect_frames_by_time(self.dirs, call.kwargs["file_prefix"])
            self.assertEqual(sorted(frames), [1., 2., 3., 4.])
        model = griddata.Model(str(self.root), "fixture")
        index = griddata.FrameIndex()
        for name in names:
            diag = diagnostics.get(name)
            base = diag.var.split(":")[0]
            for plane in ("xy", "xz"):
                with self.subTest(name=name, plane=plane):
                    prefixes = diag.candidate_prefixes("bns", plane)
                    self.assertEqual(prefixes, compare.VAR_PREFIX[base][plane])
                    expected = [2., 3., 4.] if plane == "xy" else [1., 2., 3., 4.]
                    self.assertEqual(sorted(compare.frames_by_time(self.dirs, prefixes)), expected)
                    self.assertEqual([t for t, _ in index.frames(model, diag, plane)], expected)
        self.assertIn("bhns.m1_yz", diagnostics.get("rad_E").candidate_prefixes("bhns", "yz"))

    def test_derived_fields_pair_split_and_grouped_files(self):
        for kind in ("nu_energy", "nu_avge", "nu_Ynu"):
            extra = {}
            if kind == "nu_Ynu":
                extra = dict(e_spec=("bns.m1ff_xy", "bns.rad_n_xy", "n:0"),
                             n_spec=("bns.prim_xy", "dens"))
            with self.subTest(kind=kind), \
                    patch.object(bns, "_render_avg_energy_worker_all") as worker, \
                    contextlib.redirect_stdout(io.StringIO()):
                bns.render_avg_energy_all(self.dirs, str(self.root / "figs"),
                                          derived_kind=kind, var_key=kind, **extra)
                self.assertEqual(worker.call_count, 4)
                for call in worker.call_args_list:
                    kw = call.args[0]
                    for role in ("e", "n"):
                        data = bns.BinaryData(kw[f"xy_{role}_fpath"])
                        self.assertIn(kw[f"{role}_var"], data.variable_names)
                last = worker.call_args_list[-1].args[0]
                if kind != "nu_Ynu":
                    self.assertEqual(last["xy_e_fpath"], last["xy_n_fpath"])
                    group = bns._load_frame((last["xy_e_fpath"], last["xy_n_fpath"]))
                    self.assertEqual(len(group.datasets), 1)

    def test_companion_clock_offset_and_missing_files(self):
        main = Path(self.dirs[1]) / "bin" / "bns.prim_xy.00001.bin"
        companion = main.with_name("bns.rad_n_xy.00001.bin")
        write_bin(companion, 2.04, {"n:0": 1e-5}, "xy")
        frames = bns.collect_combined_frames(
            [self.dirs[1]], "bns.prim_xy", companion_prefix=("bns.m1ff_xy", "bns.rad_n_xy"))
        self.assertEqual(frames[0][1], (str(main), str(companion)))
        self.assertEqual(bns._collect_frames_by_time(self.dirs, "missing"), {})


def write_sph(path, radii, fields, ntheta=4):
    phi, theta, radius = np.meshgrid(
        np.linspace(0, 2 * np.pi, 6, endpoint=False),
        np.linspace(.2, np.pi - .2, ntheta), radii, indexing="ij")
    points = np.stack((radius, theta, phi), axis=-1)
    npoints = radius.size
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as f:
        f.write(("# vtk DataFile Version 3.0\n"
                 "# AthenaK data at time=12 cycle=1\nBINARY\n"
                 f"DATASET STRUCTURED_GRID\nDIMENSIONS {len(radii)} {ntheta} 6\n"
                 f"POINTS {npoints} float\n").encode())
        f.write(points.astype(">f4").tobytes() + b"\n")
        f.write(f"POINT_DATA {npoints}\n".encode())
        for name, value in fields.items():
            if callable(value):
                value = value(radius, theta, phi)
            array = np.broadcast_to(value, radius.shape)
            f.write(f"SCALARS {name} float 1\nLOOKUP_TABLE default\n".encode())
            f.write(array.astype(">f4").tobytes() + b"\n")


class SphericalLayouts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        sph.set_default_radius(300)
        self.dirs = []
        for i, (grouped, radii, ntheta) in enumerate(
                ((False, [300], 4), (True, [300], 4),
                 (False, [100, 300, 500], 4), (True, [100, 300, 500], 4),
                 (True, [100, 300, 500], 5))):
            d = Path(self.tmp.name) / f"output-{i:04d}" / "sph"
            self.dirs.append(str(d))
            gauge = {"z4c_alpha": 1., "z4c_betax": 0., "z4c_betay": 0., "z4c_betaz": 0.}
            metric = {f"adm_g{axis}": float(axis in ("xx", "yy", "zz"))
                      for axis in ("xx", "yy", "zz", "xy", "xz", "yz")}
            radial = {
                "x": lambda r, t, p: .2 * np.sin(t) * np.cos(p),
                "y": lambda r, t, p: .2 * np.sin(t) * np.sin(p),
                "z": lambda r, t, p: .2 * np.cos(t),
            }
            moments = {f"{var}:{s}": value for s in range(4)
                       for var, value in (("E", 2e-9), ("N", 1e-5),
                                          *((f"F{a}", v) for a, v in radial.items()))}
            fluid = {"dens": 1e-4, "temperature": 1., "s_00": .2,
                     **{f"vel{a}": v for a, v in radial.items()}}
            groups = {"adm": metric, "mhd_w_bcc": fluid}
            if grouped:
                groups["z4c_gauge"] = gauge
                groups["m1_mom" if i == 1 else "rad_m1_moments"] = moments
            else:
                groups.update({name: {name: value} for name, value in gauge.items()})
                for field in ("E", "N", "F"):
                    groups[f"rad_m1_{field}"] = {k: v for k, v in moments.items()
                                               if k.startswith(field)}
            tag = "300.00" if len(radii) == 1 else "100.00-500.00"
            for name, fields in groups.items():
                fields = {**fields, "weights": lambda r, t, p: 4 * np.pi * r**2 / (6 * ntheta)}
                write_sph(d / f"bns.r={tag}.{name}.00001.vtk", radii, fields, ntheta)

    def test_shell_selection_and_group_resolution(self):
        styles = sph.resolve_styles(self.dirs, radius=290, var="rad_m1_F")
        self.assertEqual(len(styles), 5)
        for d, style in styles.items():
            with self.subTest(directory=d):
                path = style.path_for(d, "rad_m1_E", 1)
                data = sph.SphericalData(path, radius=290, fields=["E:0"])
                self.assertEqual(data.r, 300)
                self.assertEqual(set(data.fields), {"E:0"})
                np.testing.assert_allclose(data.fields["E:0"], 2e-9)
        with self.assertRaisesRegex(RuntimeError, "No sph VTK"):
            sph.resolve_style([self.tmp.name], var="rad_m1_F")

    def test_neutrino_results_identical_across_layouts(self):
        styles = sph.resolve_styles(self.dirs, radius=300, var="rad_m1_F")
        results = []
        fields = ("rad_m1_F", "rad_m1_E", "rad_m1_N", "adm",
                  "z4c_alpha", "z4c_betax", "z4c_betay", "z4c_betaz")
        for d in self.dirs[:4]:
            paths = tuple(styles[d].path_for(d, field, 1) for field in fields)
            results.append(neutrinos._process_snapshot((*paths, np.linspace(0, np.pi, 5))))
        self.assertTrue(np.all(np.asarray(results[0]["Lnu_E"]) > 0))
        for result in results[1:]:
            for key in results[0]:
                np.testing.assert_allclose(result[key], results[0][key], rtol=1e-12)

    def test_ejecta_results_and_changed_angular_resolution(self):
        eos = dict(nb_min=1e-10, nb_max=1., ye_min=0., ye_max=1.,
                   t_min=.1, t_max=10., h_min=1.,
                   interp_e=lambda p: np.zeros(len(p)),
                   interp_p_rho=lambda p: np.zeros(len(p)),
                   interp_s=lambda p: np.ones(len(p)))
        ejecta._init_ejecta_worker(eos, 1.67e-24, radius=300)
        styles = sph.resolve_styles(self.dirs, radius=300, var="adm")
        fields = ("mhd_w_bcc", "adm", "z4c_alpha", "z4c_betax", "z4c_betay", "z4c_betaz")
        results = []
        for d in self.dirs:
            paths = tuple(styles[d].path_for(d, field, 1) for field in fields)
            results.append(ejecta._process_snapshot_ejecta(
                (*paths, np.ones((6, 4)), np.linspace(0, np.pi, 5), 0., 100.)))
        self.assertGreater(results[0]["Mej_rate_geo"], 0)
        for result in results[1:4]:
            for key in results[0]:
                np.testing.assert_allclose(result[key], results[0][key], rtol=1e-12)
        self.assertGreater(results[-1]["Mej_rate_geo"], 0)
        self.assertEqual(results[-1]["hist2d_geo"].shape, results[0]["hist2d_geo"].shape)

    def test_extracted_shells_match_original(self):
        src = Path(self.dirs[3])
        dst = Path(self.tmp.name) / "reduced" / "sph"
        for f in sorted(src.iterdir()):
            self.assertTrue(sph_extract.extract(str(f), [100, 290], str(dst)).startswith("wrote"))
        style = sph.resolve_style([str(dst)], radius=300, var="adm")
        self.assertEqual(style.r_tag, "r=100.00-300.00")
        for f in sorted(src.iterdir()):
            reduced = style.path(str(dst), f.name.split(".")[-3], 1)
            for r in (100, 300):
                with self.subTest(file=f.name, radius=r):
                    a = sph.SphericalData(str(f), radius=r)
                    b = sph.SphericalData(reduced, radius=r)
                    self.assertEqual((b.r, b.nradii), (r, 2))
                    np.testing.assert_array_equal(b.theta, a.theta)
                    np.testing.assert_array_equal(b.phi, a.phi)
                    for name in a.fields:
                        np.testing.assert_array_equal(b.fields[name], a.fields[name])


if __name__ == "__main__":
    unittest.main()
