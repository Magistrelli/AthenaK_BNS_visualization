"""
Cut AthenaK multi-radius spherical-surface VTK dumps down to a few shells.

    python3 sph_extract.py --radii 100 200 300 400 --out-dir DEST FILE [FILE ...]

Each FILE (e.g. bns.r=10.00-500.00.adm.00184.vtk, DIMENSIONS nr ntheta nphi)
is written to DEST with only the shells nearest the requested radii, in the
same legacy-VTK layout (see sph_vtk.py), so the analysis reads it unchanged.
The output name and header carry the kept radii, e.g.
bns.r=100.00-400.00.adm.00184.vtk with nradii=4 rmin=100 rmax=400.  Values
are copied bit for bit, nothing is interpolated.  A single-radius file is
copied as it is if its radius is requested, and skipped otherwise.

Standard library only (no numpy), so it runs with any python3.  An output is
written under a temporary name, renamed once complete, and gets the input's
modification time; outputs at least as new as their input are skipped, so an
interrupted run can simply be repeated.
"""

import argparse
import os
import re
import shutil
import struct
import sys
from array import array
from multiprocessing import Pool

# bns.r=200.00.adm.12073.vtk  /  bns.r=10.00-500.00.adm.00184.vtk
_NAME_RE = re.compile(
    r'^(?P<job>.+?)\.r=(?P<rtag>[0-9]+(?:\.[0-9]+)?(?:-[0-9]+(?:\.[0-9]+)?)?)'
    r'\.(?P<rest>.+\.vtk)$')

_SIZES = {'unsigned_char': 1, 'char': 1, 'unsigned_short': 2, 'short': 2,
          'unsigned_int': 4, 'int': 4, 'float': 4,
          'unsigned_long': 8, 'long': 8, 'double': 8}
# Array codes used only to move raw items around; the byte order is irrelevant.
_WORD = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}

# Requested radii farther than this from every shell get a warning.
_R_TOL = 1e-3


class Block:
    """One section of the file: its header text and binary payload."""

    def __init__(self, head, raw=b'', itemsize=0, ncomp=0, per_point=False,
                 newline=False):
        self.head = head            # header line(s), bytes, newline included
        self.raw = raw              # binary payload
        self.itemsize = itemsize
        self.ncomp = ncomp          # values per point (per_point) or in total
        self.per_point = per_point
        self.newline = newline      # newline after the payload


def _parse(data, name):
    """Split a legacy-VTK file into Blocks; return (blocks, dims, radii)."""
    pos = 0

    def line():
        nonlocal pos
        end = data.index(b'\n', pos) + 1
        out = data[pos:end]
        pos = end
        return out

    def payload(count, itemsize):
        nonlocal pos
        raw = data[pos:pos + count * itemsize]
        if len(raw) != count * itemsize:
            raise RuntimeError(f"{name}: truncated binary block")
        pos += len(raw)
        nl = data[pos:pos + 1] == b'\n'
        pos += nl
        return raw, nl

    blocks = [Block(line() + line())]           # version banner + AthenaK meta
    if line().strip() != b'BINARY':
        raise RuntimeError(f"{name}: not a BINARY legacy VTK file")
    blocks[0].head += b'BINARY\n'
    dims = npoints = radii = None
    while pos < len(data):
        head = line()
        tok = head.split()
        if not tok:
            blocks.append(Block(head))
            continue
        key = tok[0]
        if key in (b'DATASET', b'POINT_DATA'):
            blocks.append(Block(head))
            if key == b'POINT_DATA':
                npoints = int(tok[1])
        elif key == b'DIMENSIONS':
            dims = [int(t) for t in tok[1:4]]
            blocks.append(Block(head))
        elif key == b'POINTS':
            npoints = int(tok[1])
            size = _SIZES[tok[2].decode()]
            raw, nl = payload(npoints * 3, size)
            blocks.append(Block(head, raw, size, 3, True, nl))
        elif key == b'FIELD':
            blocks.append(Block(head))
            for _ in range(int(tok[2])):
                fhead = line()
                ftok = fhead.split()
                count = int(ftok[1]) * int(ftok[2])
                size = _SIZES[ftok[3].decode()]
                raw, nl = payload(count, size)
                blocks.append(Block(fhead, raw, size, count, False, nl))
                if ftok[0] == b'RADII':
                    if ftok[3] != b'float':
                        raise RuntimeError(f"{name}: RADII is not float")
                    radii = struct.unpack(f'>{count}f', raw)
        elif key in (b'SCALARS', b'VECTORS'):
            size = _SIZES[tok[2].decode()]
            if key == b'VECTORS':
                ncomp = 3
            else:
                ncomp = int(tok[3]) if len(tok) > 3 else 1
                nxt = data[pos:data.index(b'\n', pos) + 1]
                if nxt.startswith(b'LOOKUP_TABLE'):
                    head += line()
            raw, nl = payload(npoints * ncomp, size)
            blocks.append(Block(head, raw, size, ncomp, True, nl))
        else:
            raise RuntimeError(f"{name}: unsupported VTK section {key!r}")
    if dims is None:
        raise RuntimeError(f"{name}: no DIMENSIONS")
    return blocks, dims, radii


def _select(raw, itemsize, ncomp, nr, sel):
    """Keep shells `sel` of a point array laid out as (..., nr, ncomp)."""
    code = _WORD[itemsize]
    src = array(code, raw)
    k = len(sel)
    out = array(code, bytes(len(raw) // nr * k))
    for j, ir in enumerate(sel):
        for c in range(ncomp):
            out[j * ncomp + c::k * ncomp] = src[ir * ncomp + c::nr * ncomp]
    return out.tobytes()


def _set_meta(head, **values):
    """Rewrite key=value tokens of the AthenaK metadata line."""
    banner, meta = head.split(b'\n', 1)
    for key, val in values.items():
        meta = re.sub(rb'\b' + key.encode() + rb'=\S+',
                      f'{key}={val}'.encode(), meta)
    return banner + b'\n' + meta


def extract(path, radii, out_dir):
    """Write the reduced copy of `path`; return a one-line report."""
    name = os.path.basename(path)
    m = _NAME_RE.match(name)
    if not m:
        return f"skip {name}: not an sph VTK name"
    with open(path, 'rb') as fh:
        data = fh.read()
    blocks, (nr, ntheta, nphi), file_radii = _parse(data, name)
    if file_radii is None or len(file_radii) != nr:
        pts = blocks[[b.head.startswith(b'POINTS') for b in blocks].index(True)]
        file_radii = struct.unpack(f'>{nr * 3}f', pts.raw[:nr * 3 * 4])[0::3]

    sel = sorted({min(range(nr), key=lambda i: abs(file_radii[i] - r)) for r in radii})
    warn = [r for r in radii
            if min(abs(f - r) for f in file_radii) > _R_TOL * max(1.0, r)]
    kept = [file_radii[i] for i in sel]

    if nr == 1:
        if warn:
            return f"skip {name}: r={file_radii[0]:g} not requested"
        out_name = name
    else:
        rtag = (f'r={kept[0]:.2f}-{kept[-1]:.2f}' if len(kept) > 1
                else f'r={kept[0]:.2f}')
        out_name = f"{m.group('job')}.{rtag}.{m.group('rest')}"
    out = os.path.join(out_dir, out_name)
    if os.path.exists(out) and os.path.getmtime(out) >= os.path.getmtime(path):
        return f"done {out_name}"

    tmp = out + '.tmp'
    if nr == 1:
        shutil.copyfile(path, tmp)
    else:
        npoints = len(sel) * ntheta * nphi
        with open(tmp, 'wb') as fh:
            for b in blocks:
                head, raw = b.head, b.raw
                tok = head.split() or [b'']
                if head.startswith(b'# vtk'):
                    head = _set_meta(head, nradii=len(sel),
                                     rmin=f'{kept[0]:g}', rmax=f'{kept[-1]:g}')
                elif tok[0] == b'DIMENSIONS':
                    head = f'DIMENSIONS {len(sel)} {ntheta} {nphi}\n'.encode()
                elif tok[0] in (b'POINTS', b'POINT_DATA'):
                    head = head.replace(tok[1], str(npoints).encode(), 1)
                elif tok[0] == b'RADII':
                    head = head.replace(b' ' + tok[2] + b' ',
                                        f' {len(sel)} '.encode(), 1)
                    raw = struct.pack(f'>{len(sel)}f', *kept)
                if b.per_point:
                    raw = _select(raw, b.itemsize, b.ncomp, nr, sel)
                fh.write(head + raw + (b'\n' if b.newline else b''))
    st = os.stat(path)
    os.utime(tmp, (st.st_atime, st.st_mtime))
    os.replace(tmp, out)
    msg = f"wrote {out_name}"
    if warn:
        msg += f"  (no shell at r={warn}; kept nearest {kept})"
    return msg


def _job(args):
    path, radii, out_dir = args
    try:
        return extract(path, radii, out_dir)
    except Exception as e:                      # report, keep the others going
        return f"ERROR {os.path.basename(path)}: {e}"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--radii', nargs='+', type=float, required=True,
                   help='Shell radii [M_sun] to keep (nearest shell each).')
    p.add_argument('--out-dir', required=True)
    p.add_argument('--jobs', type=int, default=1, help='Parallel processes.')
    p.add_argument('files', nargs='+')
    a = p.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    jobs = [(f, a.radii, a.out_dir) for f in a.files]
    errors = 0
    if a.jobs > 1:
        with Pool(a.jobs) as pool:
            reports = pool.imap(_job, jobs)
            for msg in reports:
                errors += msg.startswith('ERROR')
                print(msg, flush=True)
    else:
        for job in jobs:
            msg = _job(job)
            errors += msg.startswith('ERROR')
            print(msg, flush=True)
    sys.exit(1 if errors else 0)


if __name__ == '__main__':
    main()
