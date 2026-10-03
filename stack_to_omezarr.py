#!/usr/bin/env python3
"""Convert a microCT TIFF reconstruction (a folder of 2D slices) into a multiscale OME-Zarr store,
or rebuild the coarser levels of an existing OME-Zarr store from its full-resolution level.

    stack_to_omezarr.py <slice folder> <out.ome.zarr> [--levels 4] [--voxel-size UM] [--shard 128]
    stack_to_omezarr.py <in.ome.zarr> <out.ome.zarr> [--levels 4] [--shard 128]

Bruker/SkyScan NRecon folders are read as they are: the voxel size, slice range and file names
come from the *_rec.log, and the scan settings are written into the store. Any other folder of
single-page TIFFs works with --voxel-size (micrometres; slices are taken in name order).
An existing store keeps its level 0, axes, scale, origin and other attributes; only the
coarser levels are recomputed (its level 0 must be a z, y, x array).

The store is OME-Zarr 0.5 (Zarr v3): one array per level at scale{k}/<name>, 128^3 chunks, zstd,
millimetre units, origin 0 and identity LPS orientation. --shard N groups the chunks into N^3
shards (N a multiple of 128); the default, 128, writes every chunk as its own file.

Every level is the block mean of the SOURCE voxels it covers (2^k per axis), summed exactly and
rounded once (half to even) for integer data. Level k is not computed from level k-1: rounding
once keeps every coarse voxel within 0.5 gray level of the true mean, which matters for 8-bit
data, and costs no extra reading because all levels are built from the source slab in memory.

Memory: about one slab of the source, as deep as the shard, plus a coarse-level buffer (with 512^3
shards the 38 GB, 3882-slice skull peaked near 11 GB). The source is read once.
"""
import argparse
import configparser
import glob
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import tifffile
import zarr
from zarr.codecs import ZstdCodec

CHUNK = 128
ORIENTATION = {"x": "right-to-left", "y": "anterior-to-posterior", "z": "inferior-to-superior"}


def read_nrecon_log(folder):
    """(slice files, voxel size µm, acquisition attributes, name) from an NRecon *_rec.log, or None."""
    logs = glob.glob(os.path.join(folder, "*_rec.log"))
    if not logs:
        return None
    parser = configparser.ConfigParser(strict=False)
    parser.optionxform = str
    parser.read(logs[0], encoding="latin-1")
    rec, acq, system = parser["Reconstruction"], parser["Acquisition"], parser["System"]
    prefix = parser["File name convention"]["Filename Prefix"].strip()
    digits = int(parser["File name convention"].get("Filename Index Length", "8"))
    first, last = int(rec["First Section"]), int(rec["Last Section"])
    step = int(rec.get("Section to Section Step", "1"))
    extension = {"TIF": "tif", "TIFF": "tif", "BMP": "bmp", "PNG": "png"}.get(rec.get("Result File Type", "TIF").strip().upper(), "tif")
    files = [os.path.join(folder, f"{prefix}{i:0{digits}d}.{extension}") for i in range(first, last + 1, step)]
    voxel_um = float(rec["Pixel Size (um)"])
    attributes = {
        "source_log": os.path.basename(logs[0]),
        "scanner": system.get("Scanner"),
        "voxel_size_um": voxel_um,
        "source_kV": acq.get("Source Voltage (kV)", "").strip(),
        "source_uA": acq.get("Source Current (uA)", "").strip(),
        "filter": acq.get("Filter"),
        "exposure_ms": acq.get("Exposure (ms)"),
        "rotation_step_deg": acq.get("Rotation Step (deg)"),
        "study_date": acq.get("Study Date and Time"),
        "reconstruction": f'{rec.get("Reconstruction Program")} {rec.get("Program Version", "").replace("Version: ", "")}'.strip(),
        "sections": f"{first}-{last} (step {step})",
    }
    width, height = rec.get("Result Image Width (pixels)"), rec.get("Result Image Height (pixels)")
    expected = (int(height), int(width)) if width and height else None
    name = prefix.removesuffix("__rec").removesuffix("_rec").rstrip("_")
    return files, voxel_um, attributes, name, expected


def block_mean(slab, factor, dtype):
    """Mean of each factor^3 block of ``slab`` (z, y, x), trailing partial blocks dropped.
    Integers are summed exactly and rounded once, half to even."""
    z, y, x = (n // factor * factor for n in slab.shape)
    out = np.empty((z // factor, y // factor, x // factor), dtype)
    integer = np.issubdtype(dtype, np.integer)
    accumulator = np.int64 if np.issubdtype(dtype, np.signedinteger) else np.uint64 if integer else np.float64
    rows = max(factor, 32 // factor * factor)  # source slices per step: bounds the temporaries
    for z0 in range(0, z, rows):
        part = slab[z0 : min(z0 + rows, z), :y, :x]
        sums = part.reshape(part.shape[0] // factor, factor, y // factor, factor, x // factor, factor).sum(
            axis=(1, 3, 5), dtype=accumulator
        )
        mean = sums / float(factor**3)
        out[z0 // factor : z0 // factor + mean.shape[0]] = np.rint(mean) if integer else mean
    return out


class StackSource:
    """Slabs of a TIFF slice folder (NRecon log read when present)."""

    def __init__(self, folder, voxel_um=None, threads=8):
        nrecon = read_nrecon_log(folder)
        if nrecon:
            files, log_voxel_um, attributes, name, expected = nrecon
            voxel_um = voxel_um or log_voxel_um
        else:
            if voxel_um is None:
                sys.exit("No *_rec.log in the folder: give --voxel-size in micrometres")
            files = sorted(f for f in glob.glob(os.path.join(folder, "*")) if f.lower().endswith((".tif", ".tiff")))
            attributes, name, expected = {"voxel_size_um": voxel_um}, os.path.basename(os.path.normpath(folder)), None
        missing = [f for f in files if not os.path.exists(f)]
        if missing:
            sys.exit(f"{len(missing)} slices missing, e.g. {missing[0]}")
        page = tifffile.TiffFile(files[0]).pages[0]
        if expected and expected != tuple(page.shape):
            sys.exit(f"slices are {page.shape[0]}x{page.shape[1]}, the log says {expected[0]}x{expected[1]}")
        self.files, self.name, self.dtype = files, name, page.dtype
        self.shape = (len(files), *page.shape)
        self.spacing = [voxel_um / 1000.0] * 3  # mm, z y x
        self.origin = [0.0] * 3
        self.axes = [{"name": a, "type": "space", "unit": "millimeter",
                      "orientation": {"type": "anatomical", "value": ORIENTATION[a]}} for a in "zyx"]
        self.attributes = {"acquisition": attributes}
        self.pool = ThreadPoolExecutor(threads)
        self.description = f"{voxel_um} µm isotropic"

    def slab(self, z0, z1):
        slab = np.empty((z1 - z0, *self.shape[1:]), dtype=self.dtype)

        def load(i):
            slab[i] = tifffile.imread(self.files[z0 + i])

        list(self.pool.map(load, range(z1 - z0)))
        return slab

    def reference(self, z):
        return tifffile.imread(self.files[z])


class StoreSource:
    """Level 0 of an existing OME-Zarr store, with its axes, transforms and other attributes."""

    def __init__(self, path):
        group = zarr.open_group(path, mode="r")
        attrs = dict(group.attrs)
        multiscales = attrs["ome"]["multiscales"][0]
        dataset = multiscales["datasets"][0]
        self.array = group[dataset["path"]]
        names = [a["name"] for a in multiscales["axes"]]
        if names != ["z", "y", "x"]:
            sys.exit(f"level 0 axes are {names}; only z, y, x stores are handled")
        transforms = {t["type"]: t for t in dataset["coordinateTransformations"]}
        self.spacing = [float(v) for v in transforms["scale"]["scale"]]
        self.origin = [float(v) for v in transforms.get("translation", {"translation": [0.0] * 3})["translation"]]
        self.axes = multiscales["axes"]
        self.name = multiscales.get("name") or dataset["path"].split("/")[-1]
        self.shape, self.dtype = tuple(self.array.shape), self.array.dtype
        self.attributes = {k: v for k, v in attrs.items() if k != "ome"}
        self.description = f"level 0 of {os.path.basename(os.path.normpath(path))}, spacing {self.spacing} mm"

    def slab(self, z0, z1):
        return self.array[z0:z1]

    def reference(self, z):
        return self.array[z]


def convert(source, out, levels=4, shard=CHUNK, log=print):
    """Write ``source`` (StackSource or StoreSource) to ``out`` with ``levels`` levels, chunks grouped in
    ``shard``^3 shards (``shard`` == CHUNK: no shards, one file per chunk)."""
    if shard < CHUNK or shard % CHUNK:
        raise ValueError(f"shard size {shard} is not a multiple of the {CHUNK}-voxel chunk")
    name, dtype = source.name, source.dtype
    log(f"{name}: {source.shape[2]}x{source.shape[1]}x{source.shape[0]} {dtype}, {source.description}, {levels} levels")
    shapes = [tuple(n // 2**k for n in source.shape) for k in range(levels)]
    root = zarr.open_group(out, mode="w", zarr_format=3)
    arrays = []
    for k, shape in enumerate(shapes):
        group = root.create_group(f"scale{k}")
        group.attrs["_ARRAY_DIMENSIONS"] = ["z", "y", "x"]
        layout = {"chunks": (CHUNK,) * 3, "shards": (shard,) * 3} if shard > CHUNK else {"chunks": (CHUNK,) * 3}
        arrays.append(
            group.create_array(
                name, shape=shape, dtype=dtype, compressors=ZstdCodec(level=3), fill_value=0,
                dimension_names=("z", "y", "x"), config={"write_empty_chunks": False}, **layout,
            )
        )
    datasets = []
    for k in range(levels):
        factor = 2**k
        datasets.append({
            "path": f"scale{k}/{name}",
            "coordinateTransformations": [
                {"type": "scale", "scale": [s * factor for s in source.spacing]},
                # a level-k voxel is the mean of a factor^3 block: its centre sits (factor-1)/2 voxels in
                {"type": "translation", "translation": [o + s * (factor - 1) / 2.0 for o, s in zip(source.origin, source.spacing)]},
            ],
        })
    root.attrs["ome"] = {
        "version": "0.5",
        "multiscales": [{
            "name": name,
            "axes": source.axes,
            "datasets": datasets,
            "type": "local_mean",
            "metadata": {
                "description": "Each level is the mean of the 2^k x 2^k x 2^k level-0 voxels it covers, rounded once (half to even)",
                "method": "stack_to_omezarr.block_mean",
            },
        }],
    }
    for key, value in source.attributes.items():
        root.attrs[key] = value

    pending = [[] for _ in range(levels)]  # coarse slices waiting to fill a row of shards
    written = [0] * levels
    # whole shards per write, and a multiple of every 2^k so blocks never straddle slabs
    rows = math.lcm(shard, 2 ** (levels - 1))

    def flush(k, final=False):
        data = np.concatenate(pending[k]) if pending[k] else None
        while data is not None and (len(data) >= rows or (final and len(data))):
            n = min(rows, len(data))
            arrays[k][written[k] : written[k] + n] = data[:n]
            written[k] += n
            data = data[n:] if len(data) > n else None
        pending[k] = [data] if data is not None else []

    started = time.time()
    depth = source.shape[0]
    for z0 in range(0, depth, rows):
        z1 = min(z0 + rows, depth)
        slab = np.asarray(source.slab(z0, z1))
        arrays[0][z0:z1] = slab
        written[0] += z1 - z0
        for k in range(1, levels):
            pending[k].append(block_mean(slab, 2**k, dtype))
            flush(k)
        del slab
        log(f"  slices {z0}-{z1 - 1} done, {time.time() - started:.0f} s")
    for k in range(1, levels):
        flush(k, final=True)
    assert written == [s[0] for s in shapes], (written, shapes)
    zarr.consolidate_metadata(out)
    log(f"levels {shapes} written in {time.time() - started:.0f} s")


def verify(out, source, levels, log=print, samples=6, seed=0):
    """Level 0 against the source; every coarse level against the exact block mean of level 0."""
    rng = np.random.default_rng(seed)
    name, depth = source.name, source.shape[0]
    group = zarr.open_group(out, mode="r")
    level0 = group[f"scale0/{name}"]
    mismatched = 0
    for z in rng.choice(depth, min(samples, depth), replace=False):
        mismatched += int((level0[z] != np.asarray(source.reference(z))).sum())
    log(f"level 0 vs source, {min(samples, depth)} random slices: {mismatched} voxels differ")
    ok = mismatched == 0
    for k in range(1, levels):
        factor, array = 2**k, group[f"scale{k}/{name}"]
        worst = 0.0
        for _ in range(3):
            n = min(8, *array.shape)
            cz, cy, cx = (int(rng.integers(0, s - n + 1)) for s in array.shape)
            source = level0[cz * factor : (cz + n) * factor, cy * factor : (cy + n) * factor, cx * factor : (cx + n) * factor]
            mean = source.astype(np.float64).reshape(n, factor, n, factor, n, factor).mean(axis=(1, 3, 5))
            stored = array[cz : cz + n, cy : cy + n, cx : cx + n].astype(np.float64)
            worst = max(worst, float(np.abs(stored - mean).max()))
        log(f"level {k} {array.shape}: max |stored - true {factor}^3 mean| = {worst:.3f} gray levels")
        ok &= worst <= 0.5 + 1e-9
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("folder", help="folder of 2D TIFF slices (an NRecon *_Rec folder is read with its log), or an .ome.zarr store")
    parser.add_argument("out", help="output .ome.zarr directory (replaced if it exists)")
    parser.add_argument("--levels", type=int, default=4, help="resolution levels including full resolution (default 4)")
    parser.add_argument("--voxel-size", type=float, help="voxel size in micrometres (overrides the log)")
    parser.add_argument("--threads", type=int, default=8, help="threads reading slices (default 8)")
    parser.add_argument("--shard", type=int, default=CHUNK,
                        help=f"shard size in voxels per side, a multiple of {CHUNK} (default {CHUNK}: no shards, one file per chunk)")
    parser.add_argument("--no-verify", action="store_true", help="skip the check against the source")
    args = parser.parse_args()
    if args.shard < CHUNK or args.shard % CHUNK:
        parser.error(f"--shard must be a multiple of {CHUNK}")
    if os.path.samefile(args.folder, args.out) if os.path.exists(args.out) else False:
        sys.exit("write to a new store: the source is read while the output is written")
    if os.path.exists(os.path.join(args.folder, "zarr.json")):
        source = StoreSource(args.folder)
    else:
        source = StackSource(args.folder, args.voxel_size, args.threads)
    convert(source, args.out, args.levels, shard=args.shard)
    if not args.no_verify and not verify(args.out, source, args.levels):
        sys.exit("VERIFY FAILED")
    print("DONE")


if __name__ == "__main__":
    main()
