#pragma once
#include "slabcc_math.hpp"

struct slabcc_cell {
	rowvec3 vec_lengths = { 0, 0, 0 };	// length of the basis vectors in Bohr
	mat33 vectors = zeros(3, 3);		// supercell vectors in Bohr
	urowvec3 grid = { 1, 1, 1 };		// grid size
	uword normal_direction = 0;
	double voxel_vol = 0;				// voxel volume of the grid in Bohr^3

	//Sets the global struct slabcc_cell parameters from the new cell vectors (in Bohr) and the grid density
	void update(const mat33& new_vectors, const urowvec3& new_grid);
};
