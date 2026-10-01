"""Synthetic checks: python -m pytest test_stack_to_omezarr.py (or run this file directly)."""
import os
import tempfile

import numpy as np
import tifffile
import zarr

import stack_to_omezarr as conv

LOG = """[System]
Scanner=SkyScan1272
[Acquisition]
Source Voltage (kV)=  70
Study Date and Time=28 Aug 2023  22h:25m:41s
[Reconstruction]
Reconstruction Program=NRecon
First Section=5
Last Section={last}
Section to Section Step=1
Result File Type=TIF
Result Image Width (pixels)={width}
Result Image Height (pixels)={height}
Pixel Size (um)=6.12669
[File name convention]
Filename Index Length=8
Filename Prefix=synthetic__rec
"""


def make_stack(folder, dtype, shape, seed=0):
    rng = np.random.default_rng(seed)
    depth, height, width = shape
    top = np.iinfo(dtype).max
    volume = rng.integers(0, top, size=shape, endpoint=True, dtype=dtype)
    for i in range(depth):
        tifffile.imwrite(os.path.join(folder, f"synthetic__rec{i + 5:08d}.tif"), volume[i])
    with open(os.path.join(folder, "synthetic__rec.log"), "w") as log:
        log.write(LOG.format(last=depth + 4, width=width, height=height))
    return volume


def check(dtype, shape, levels=4):
    with tempfile.TemporaryDirectory() as tmp:
        source = os.path.join(tmp, "synthetic_Rec")
        os.mkdir(source)
        volume = make_stack(source, dtype, shape)
        out = os.path.join(tmp, "synthetic.ome.zarr")
        stack = conv.StackSource(source)
        conv.convert(stack, out, levels=levels, log=lambda *a: None)
        name = stack.name
        assert name == "synthetic"
        group = zarr.open_group(out, mode="r")
        np.testing.assert_array_equal(group[f"scale0/{name}"][:], volume)
        for k in range(1, levels):
            f = 2**k
            z, y, x = (n // f for n in volume.shape)
            exact = volume[: z * f, : y * f, : x * f].astype(np.float64).reshape(z, f, y, f, x, f).mean(axis=(1, 3, 5))
            stored = group[f"scale{k}/{name}"][:]
            assert stored.shape == (z, y, x)
            np.testing.assert_array_equal(stored, np.rint(exact).astype(dtype))  # rounded once, from the source
        ms = group.attrs["ome"]["multiscales"][0]
        assert ms["datasets"][2]["coordinateTransformations"][0]["scale"] == [0.00612669 * 4] * 3
        assert ms["datasets"][2]["coordinateTransformations"][1]["translation"] == [0.00612669 * 1.5] * 3
        assert group.attrs["acquisition"]["voxel_size_um"] == 6.12669
        assert conv.verify(out, stack, levels, log=lambda *a: None)

        # Rebuilding from the store itself (plain chunks this time) gives the same levels and keeps its metadata.
        again = os.path.join(tmp, "again.ome.zarr")
        store = conv.StoreSource(out)
        conv.convert(store, again, levels=levels, sharded=False, log=lambda *a: None)
        rebuilt = zarr.open_group(again, mode="r")
        for k in range(levels):
            np.testing.assert_array_equal(rebuilt[f"scale{k}/{name}"][:], group[f"scale{k}/{name}"][:])
        assert rebuilt.attrs["ome"]["multiscales"][0]["datasets"] == ms["datasets"]
        assert rebuilt.attrs["acquisition"] == group.attrs["acquisition"]
        assert "sharding_indexed" not in str(rebuilt[f"scale0/{name}"].metadata.codecs)
        assert conv.verify(again, store, levels, log=lambda *a: None)


def test_uint8_odd_shape():
    check(np.uint8, (523, 37, 41))  # crosses a 512-slice slab; odd sizes drop trailing partial blocks


def test_uint16():
    check(np.uint16, (40, 66, 70))


if __name__ == "__main__":
    test_uint8_odd_shape()
    test_uint16()
    print("all passed")
