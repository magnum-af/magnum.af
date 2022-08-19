#!/usr/bin/python3
# Example demonstrating usage of regions

import arrayfire as af

# helper functions


def num_cells_in_region(regions, select_region):
    """Number of cells which match the select_region index."""
    one_if_index_matches_else_zero = (regions == select_region)
    num_of_matching_cells = af.sum(one_if_index_matches_else_zero)
    return num_of_matching_cells


def mean_in_region(m, regions, select_region):
    """Returns the spacially averaged values of a vector field evaluated in the selected region only."""
    mean = af.sum(af.sum(af.sum(m * af.tile(regions == select_region, 1, 1, 1, 3),
                  dim=0), dim=1), dim=2)/num_cells_in_region(regions, select_region)
    return [mean[:, :, :, 0].scalar(), mean[:, :, :, 1].scalar(), mean[:, :, :, 2].scalar()]


def field_in_region(field, regions, select_region):
    """Returns the field in the selected region, all other cells are set to zero. This can be used e.g. to write the region to a file."""
    return field * af.tile(regions == select_region, 1, 1, 1, 3)


# setting size
nx, ny, nz = 2, 3, 3

# regions are defined per cell. Size is [nx, ny, nz, 1].
regions = af.constant(0., nx, ny, nz, 1, dtype=af.Dtype.u16)
regions[:, :, 0] = 1
regions[:, :, 1] = 2
regions[:, :, 2] = 3
print(regions)


# setting m
m = af.constant(0., nx, ny, nz, 3, dtype=af.Dtype.f64)
m[:, :, 0, 0] = 1.  # i.e. mx = 1 in layer 0
m[:, :, 1, 1] = 1.  # i.e. my = 1 in layer 1
m[:, :, 2, 2] = 1.  # i.e. mz = 1 in layer 2
print(m)

# get the magnetization in one region only, else is zero
print(field_in_region(m, regions, select_region=1))
print(field_in_region(m, regions, select_region=2))
print(field_in_region(m, regions, select_region=3))

# evaluating mean in regions only
print(mean_in_region(m, regions, select_region=1))  # should give [1, 0, 0]
print(mean_in_region(m, regions, select_region=2))  # should give [0, 1, 0]
print(mean_in_region(m, regions, select_region=3))  # should give [0, 0, 1]

# sanity checks
print(num_cells_in_region(regions, select_region=1))
print(num_cells_in_region(regions, select_region=2))
print(num_cells_in_region(regions, select_region=3))
