// Copyright (c) 2018, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#pragma once
#include "slabcc_math.hpp"

using namespace std;
extern const double ang_to_bohr;

// These functions are processing the VASP files:
// read/write POSCAR, CHGCAR, LOCPOT files
// convert direct <> relative coordinates
// swap_axes of the model and all its loaded data
// calculate the planar average of the 3D grid-based data


//missing the spin polarized containers
//NOT SUITABLE FOR GENERAL PURPOSE APPLICATIONS!
struct supercell {
	struct atom {

		//the order in which the atoms are declated in the POSCAR
		//needed to keep the atom groups with the same name seperated
		vector<int> definition_order;

		//atomic types as declared in the POSCAR
		vector<string> type;

		// all atomic positions
		mat position;

		//relaxation constraints for each atom defined as "T" or "F" in the POSCAR
		//if selective_dynamics is false, these must be ignored
		vector <vector<bool>> constrains;
	};

	//VASP POSCAR label
	string label;

	// VASP POSCAR file scaling factor
	double scaling = 1;

	//VASP cell vectors in the POSCAR
	mat33 cell_vectors;

	bool selective_dynamics = true;

	//cartesian: if definition line in the POSCAR starts with "c" or "k"
	//direct:	if definition line in the POSCAR starts with "d"
	string coordination_system;

	atom atoms;

	//total number of atoms
	int atoms_number = 0;

	cube charge;			//VASP first dataset (spin 1+2) in CHGCAR
	cube potential;			//VASP first dataset (spin 1+2) in LOCPOT
};


//Returns the direct (relative) position based on the provided cartesians and the supercell size*scaling
rowvec3 direct_cord(const supercell& structure, const rowvec3& cartesians);

//normalizes the coordinates in direct/relative system to [0 1]
void normalize_positions(supercell& structure);

//generates a supercell from the POSCAR file
supercell read_POSCAR(const string& file_name);

//reads grid data from CHGCAR/LOCPOT files and return ONLY the first data set (spin1+2)
//NOT SUITABLE FOR GENERAL PURPOSE APPLICATIONS!
cube read_CHGPOT(const string& file_name);

//shifts the whole supercell (positions, charge, potential) by pos_shift as relative shift vector [0 1]
void shift_structure(supercell& structure, const rowvec3& pos_shift);

//type: "CHGCAR", "LOCPOT"
void write_CHGPOT(const string& type, const string& file_name, const supercell& structure);

//writes POSCAR of a supercell to a file
void write_POSCAR(const supercell& structure, const string& file_name);

//Planar average of grid-based data in the desired direction from a supecell
//type: "CHGCAR", "POTCAR"
//direction: 0,1,2 > x,y,z
//correct: apply the normalization
vector<double> planar_average(const string& type, const uword& direction, const bool& correct, const supercell& structure);

void write_planar_avg(const supercell& structure, const string& id);

void write_planar_avg(const cx_cube& potential_data, const cx_cube& charge_data, const string& id);
