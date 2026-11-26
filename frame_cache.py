"""
Bookkeeping of a rendered frame series, so a rerun only renders what changed.

A series lives in one folder (e.g. Figs/dens_xy/).  Its frames are PNG files in
the folder, or -- once make_movies.sh / movie_tools.sh has encoded the movie --
packed into frames.tar there.  frames.json records, for every frame, the input
files it was rendered from (name, size, modification time) and the plot
settings of the whole series.  A frame counts as done if it exists (on disk or
in frames.tar) and its inputs are unchanged; if the settings changed, the
whole series starts over.  A folder from before this bookkeeping (no
frames.json) keeps its existing frames, as the old skip-existing did.

Standard library only.
"""

import json
import os
import tarfile

TAR = "frames.tar"
MANIFEST = "frames.json"


def stamp(paths):
    """(name, size, mtime) of every input file; renamed copies with the same
    content and times (e.g. an archived run) keep the same stamp."""
    flat = []
    for p in paths:
        if isinstance(p, (list, tuple)):
            flat.extend(p)
        elif p:
            flat.append(p)
    out = []
    for p in flat:
        try:
            st = os.stat(p)
            out.append([os.path.basename(p), st.st_size, int(st.st_mtime)])
        except OSError:
            out.append([os.path.basename(p), None, None])
    return out


class FrameSeries:
    def __init__(self, folder, settings):
        self.folder = folder
        os.makedirs(folder, exist_ok=True)
        # json round trip, so tuples and lists compare equal to the stored ones
        self.settings = json.loads(json.dumps(settings, default=repr))
        self.manifest = os.path.join(folder, MANIFEST)
        self.tar = os.path.join(folder, TAR)
        self.legacy = not os.path.exists(self.manifest)
        self.frames = {}
        if not self.legacy:
            with open(self.manifest) as fh:
                saved = json.load(fh)
            if saved.get("settings") == self.settings:
                self.frames = saved.get("frames", {})
            else:
                self.restart()
        self.packed = set()
        if os.path.exists(self.tar):
            with tarfile.open(self.tar) as tf:
                self.packed = set(tf.getnames())

    def restart(self):
        """Settings changed: drop every frame of the series."""
        for name in os.listdir(self.folder):
            if name.endswith(".png") or name == TAR:
                os.remove(os.path.join(self.folder, name))
        self.frames = {}

    def done(self, name, sources):
        """True if frame `name` exists and was rendered from these sources."""
        if not (os.path.exists(os.path.join(self.folder, name)) or name in self.packed):
            return False
        if self.legacy and name not in self.frames:
            self.frames[name] = stamp(sources)      # trust frames from before
            return True
        return self.frames.get(name) == stamp(sources)

    def record(self, name, sources):
        if os.path.exists(os.path.join(self.folder, name)):
            self.frames[name] = stamp(sources)

    def save(self):
        tmp = self.manifest + ".tmp"
        with open(tmp, "w") as fh:
            json.dump({"settings": self.settings, "frames": self.frames}, fh, indent=0)
        os.replace(tmp, self.manifest)
