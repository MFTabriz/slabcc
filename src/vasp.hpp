// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#pragma once
#include "slabcc_math.hpp"
// TODO: CONSTEXPER
extern const double ang_to_bohr;

// These functions are processing the VASP files:
// read/write POSCAR, CHGCAR, LOCPOT files
// convert direct <> relative coordinates
// calculate the planar average of the 3D grid-based data

// missing the spin polarized containers
// NOT SUITABLE FOR GENERAL PURPOSE APPLICATIONS!
struct supercell {
  struct atom {

    // the order in which the atoms are declared in the POSCAR
    // needed to keep the atom groups with the same name separated
    std::vector<int> definition_order;

    // atomic types as declared in the POSCAR
    std::vector<std::string> type;

    // all atomic positions
    arma::mat position;

    // relaxation constraints for each atom defined as "T" or "F" in the POSCAR
    // if selective_dynamics is false, these must be ignored
    std::vector<std::vector<bool>> constrains;
  };

  // VASP POSCAR label
  std::string label;

  // VASP POSCAR file scaling factor
  double scaling = 1;

  // VASP cell vectors in the POSCAR
  arma::mat33 cell_vectors;

  bool selective_dynamics = true;

  // cartesian: if definition line in the POSCAR starts with "c" or "k"
  // direct: if definition line in the POSCAR starts with "d"
  std::string coordination_system;

  atom atoms;

  // total number of atoms
  arma::uword atoms_number = 0;

  // total charge distribution, negative for presence of electron
  // (spin 1+2 dataset of VASP CHGCAR, * -1 / volume (Bohr^3))
  arma::cube charge;

  arma::cube potential; // total potential (VASP LOCPOT * -1)

  // generates a supercell and loads its data the POSCAR file
  explicit supercell(const std::string &file_name);

  // shifts the whole supercell (positions, charge, potential) by pos_shift as
  // relative shift vector [0 1]
  void shift(const arma::rowvec3 &pos_shift);

  // Returns the direct (relative) position based on the provided Cartesians and
  // the supercell size*scaling
  arma::rowvec3 direct_cord(const arma::rowvec3 &cartesians);

  void write_CHGCAR(const std::string &file_name) const;
  void write_LOCPOT(const std::string &file_name) const;

private:
  // type: "CHGCAR", "LOCPOT"
  void write_CHGPOT(const std::string &type,
                    const std::string &file_name) const;

  // writes POSCAR of a supercell to a file
  void write_POSCAR(const std::string &file_name) const;

  // normalizes the coordinates in direct/relative system to [0 1]
  void normalize_positions();
};

// reads grid data from CHGCAR/LOCPOT files and return ONLY the first data set
// (spin 1+2) NOT SUITABLE FOR GENERAL PURPOSE APPLICATIONS!
arma::cube read_VASP_grid_data(const std::string &file_name);

// Write planar average of potential and charge to files
// coordinate_vectors (Bohr)
void write_planar_avg(const arma::cube &potential_data,
                      const arma::cube &charge_data, const std::string &id,
                      const arma::rowvec3 &coordinate_vectors,
                      const int direction = -1);

// check conditions and consistency of the supercell grid sizes and the shape
void check_slabcc_compatiblity(const supercell &Neutral_supercell,
                               const supercell &Charged_supercell);