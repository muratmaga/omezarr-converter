# omezarr-converter

Converts a microCT reconstruction (a folder of 2D TIFF slices) into a multiscale
[OME-Zarr](https://ngff.openmicroscopy.org/) store that can be streamed from object storage. It
can also rebuild the coarser levels of an existing OME-Zarr store from its full-resolution level.

Bruker/SkyScan NRecon output folders are read as they are. The voxel size, slice range and file
names come from the `*_rec.log`, and the scan settings are written into the store.

## Usage

```bash
pip install -r requirements.txt

# NRecon reconstruction folder (voxel size from the *_rec.log)
python stack_to_omezarr.py Scan_0006_Rec Scan_0006.ome.zarr

# any other folder of single-page TIFFs (slices taken in name order)
python stack_to_omezarr.py slices/ specimen.ome.zarr --voxel-size 12.5

# rebuild the levels of an existing store into a new one (level 0, axes, origin kept)
python stack_to_omezarr.py old.ome.zarr new.ome.zarr

# group the chunks into 384³ shards (fewer files, for very large volumes or object storage)
python stack_to_omezarr.py Scan_0006_Rec Scan_0006.ome.zarr --shard 384
```

| Option | Default | |
|---|---|---|
| `--levels` | 4 | resolution levels, full resolution included |
| `--voxel-size` | from the log | micrometres, isotropic |
| `--shard` | 128 | shard size in voxels per side, a multiple of 128. At 128 every 128³ chunk is its own file (no shards) |
| `--threads` | 8 | threads reading slices |
| `--no-verify` | off | skip the check against the source |

After writing, the store is checked against its source. Level 0 must match the source slices
exactly, and every coarse level must be within 0.5 gray level of the exact block mean. A failed
check exits with `VERIFY FAILED`.

## Output

- OME-Zarr 0.5 (Zarr v3). Each level is one array at `scale{k}/<name>`, in z, y, x order.
- 128³ chunks, compressed with zstd. Chunks that are all zero are not written.
- With `--shard N` the chunks are grouped into N³ shards: one file per shard instead of one per
  chunk. Shards matter when the store goes to object storage (each file read is one request) or
  sits on a network drive or a synced folder. On a local disk plain chunks are the simpler choice.
- Millimetre units. The origin is 0 and the orientation is identity LPS, because a slice folder carries neither.
- Scan settings are stored in the group attribute `acquisition`.

## How the coarse levels are computed

Each voxel of level *k* is the mean of the 2^k × 2^k × 2^k level-0 voxels it covers. It is
computed from level 0, not from level *k*−1. Integer sums are exact and rounded once (half to
even). Every coarse voxel is therefore within 0.5 gray level of the true mean, which matters for
8-bit data. Building each level from the one above it would round at every step, reaching about
±1 gray level by level 3.

This costs no extra reading. The source is read once, in slabs as deep as the shard (128 slices
by default), and every level is built from the slab in memory. Peak memory is about one slab
plus a coarse-level buffer, so larger shards need more memory. With 512³ shards a
1908 × 2596 × 3882 16-bit volume (38 GB raw) peaked near 11 GB and converted in about 2 minutes.

Each level's voxel centre sits (2^k − 1)/2 level-0 voxels in from the origin, and the
`translation` of each level records this. Every level therefore lines up in physical space.

## Tests

```bash
python -m pytest test_stack_to_omezarr.py   # or: python test_stack_to_omezarr.py
```

The tests build synthetic 8-bit and 16-bit stacks, including odd sizes and a stack that crosses
several slabs, with and without shards. They check every coarse level against the exactly rounded
mean, and check that rebuilding from a store in the other layout gives identical levels and keeps
its metadata.
