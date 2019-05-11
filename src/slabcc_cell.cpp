// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "slabcc_cell.hpp"

void slabcc_cell::update(const mat33& new_vectors, const urowvec3& new_grid) {
	vectors = new_vectors;
	grid = new_grid;
	for (uword i = 0; i < 3; ++i) {
		vec_lengths(i) = norm(vectors.col(i));
	}
	voxel_vol = prod(vec_lengths / grid);
}
