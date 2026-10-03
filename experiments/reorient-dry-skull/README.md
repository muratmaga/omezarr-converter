# Reorient + crop of the dry skull (inputs)

These are the inputs for the experiment tracked in issue #1. They were made in 3D Slicer's Crop Volume module ("Reorient" section) on level 1 of the dry-skull store and are meant to be applied to level 0.

- `Reorient_transform.h5`: ITK `AffineTransform_double_3_3`, LPS. This is Slicer's "from parent" transform, so the rotation applied to the volume is its inverse. It is a rigid rotation of 138°, with no scaling.
- `Crop_Extra.mrk.json`: a Slicer Markups ROI, aligned to the world axes. Center (LPS) 6.467, 7.558, 11.654 mm; size 14.954 × 26.087 × 11.811 mm. The box was placed by hand.

## Data

`https://js2.jetstream-cloud.org:8001/MorphoDepot-volumes/large-sample-data/dry_skull_4K_12.5micron.ome.zarr`. Level 0 is 1908 × 2596 × 3882 (x, y, z) uint16, 6.12669 µm (the name says 12.5 µm, but the scan is 6.127 µm). Every level shares one physical frame, so the transform and box apply to any level unchanged.

## Replicating it in Slicer

1. Load level 1 (`scale1`) of the store.
2. Load both files and put the volume under the transform.
3. In Crop Volume, harden the transform: Reorient → Apply. This only changes the volume's orientation matrix; no voxels are touched.
4. Crop with the ROI using interpolated cropping (spacing scale 1). This step resamples every voxel.

Resulting size at level 1: 1220 × 2129 × 963 voxels.
