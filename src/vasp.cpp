// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "vasp.hpp"

void supercell::write_POSCAR(const std::string &file_name) const {
  std::ofstream out_file;
  out_file.open(file_name);
  out_file << label << '\n';
  out_file << "   " << std::fixed << std::showpos << std::setprecision(16)
           << scaling << '\n';
  cell_vectors.t().raw_print(out_file);
  std::vector<int> defs = atoms.definition_order;
  const auto end =
      std::unique(defs.begin(), defs.end(),
                  [](const int &l, const int &r) noexcept { return l == r; });
  defs.erase(end, defs.end());
  for (const auto &i : defs) {
    const auto it = std::find(atoms.definition_order.begin(),
                              atoms.definition_order.end(), i);
    const auto pos = std::distance(atoms.definition_order.begin(), it);
    out_file << " " << atoms.type.at(pos);
  }
  out_file << '\n' << std::noshowpos;

  arma::uword counter = 0;
  arma::uword counted = 0;
  for (size_t k = 0; k < defs.size(); ++k) {
    while (atoms.definition_order.at(counter) == defs.at(k)) {
      ++counter;
      if (counter == atoms_number)
        break;
    }
    out_file << " " << counter - counted;
    counted = counter;
  }
  out_file << '\n';

  if (selective_dynamics) {
    out_file << "Selective dynamics\n";
  }
  out_file << coordination_system << '\n';

  for (arma::uword number = 0; number < atoms_number; ++number) {
    out_file << " " << atoms.position.row(number);
    if (selective_dynamics) {
      for (auto qq2 = 0; qq2 < 3; ++qq2) {
        if (atoms.constrains.at(number).at(qq2)) {
          out_file << " F";
        } else {
          out_file << " T";
        }
      }
    }
    out_file << '\n';
  }

  out_file.close();
}

arma::rowvec3 supercell::direct_cord(const arma::rowvec3 &cartesians) {
  const arma::mat33 scaled_cell_vectors = cell_vectors / scaling;
  const arma::rowvec3 direct_coords =
      arma::solve(scaled_cell_vectors.t(), cartesians.t()).t();
  return direct_coords;
}

void supercell::normalize_positions() {
  auto log = spdlog::get("loggers");
  if (coordination_system == "direct") {
    atoms.position = fmod_p(atoms.position, 1);
  } else {
    // It can be converted to direct coordinates, normalized and converted back.
    log->warn("Position normalizion in the Cartesian coordinates has not been "
              "implemented yet!");
  }
}

void supercell::shift(const arma::rowvec3 &relative_shift) {

  arma::rowvec3 pos_shift = relative_shift;

  // if our atomic positions are in cartesian coordinates,
  // find the equivalent cartesian vector
  if (coordination_system == "cartesian") {
    pos_shift = direct_cord(pos_shift);
  }

  atoms.position.each_row() += pos_shift;

  if (coordination_system == "direct") {
    normalize_positions();
  }

  charge = ::shift(charge, relative_shift);
  potential = ::shift(potential, relative_shift);
}

supercell::supercell(const std::string &file_name) {
  auto log = spdlog::get("loggers");
  std::ifstream infile;
  std::string temp_line;
  infile.open(file_name);
  if (!infile) {
    log->critical("Could not open the " + file_name);
  }
  // TODO: if file is unreadable, message!
  std::getline(infile, label);
  std::getline(infile, temp_line);
  scaling = stod(temp_line);
  infile >> cell_vectors;

  std::getline(infile, temp_line);
  std::getline(infile, temp_line);

  std::string buf;
  std::stringstream ss(temp_line);
  std::vector<std::string> temp_types;
  std::vector<int> temp_numbers;

  while (ss >> buf) {
    temp_types.push_back(buf);
  }
  std::getline(infile, temp_line);
  std::stringstream ss2(temp_line);
  while (ss2 >> buf) {
    temp_numbers.push_back(stoi(buf));
  }

  int atom_counter = 0;
  for (size_t i = 0; i < temp_types.size(); ++i) {
    for (auto a = 0; a < temp_numbers.at(i); ++a) {
      atoms.type.push_back(temp_types.at(i));
      atoms.definition_order.push_back(i);
      ++atom_counter;
    }
  }
  atoms_number = atom_counter;
  atoms.position.set_size(atom_counter, 3);
  atoms.constrains.resize(atom_counter, std::vector<bool>(3, 0));
  std::getline(infile, temp_line);
  if (tolower(temp_line.at(0)) == 's') {
    selective_dynamics = true;
    std::getline(infile, temp_line);
  } else {
    selective_dynamics = false;
  }

  if ((tolower(temp_line.at(0)) == 'k') || (tolower(temp_line.at(0)) == 'c')) {
    coordination_system = "cartesian";
  } else if (tolower(temp_line.at(0)) == 'd') {
    coordination_system = "direct"; // in VASP it means relative!!
  } else {
    log->critical("Unknown coordination system type in " + file_name);
  }

  for (arma::uword i = 0; i < atoms_number; ++i) {
    std::string temp_constrain;
    infile >> atoms.position.row(i);
    if (selective_dynamics) {
      for (auto qq2 = 0; qq2 < 3; ++qq2) {
        infile >> temp_constrain;
        atoms.constrains.at(i).at(qq2) =
            (tolower(temp_constrain.at(0)) == 't') ? false : true;
      }
    } else {
      infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
  }

  normalize_positions();
}

arma::cube read_VASP_grid_data(const std::string &file_name) {
  auto log = spdlog::get("loggers");
  std::ifstream infile;
  std::string temp_line;
  infile.open(file_name);
  if (!infile) {
    log->critical("File not found: " + file_name);
    return {};
  }
  const supercell structure(file_name);
  for (arma::uword currLineNumber = 0;
       currLineNumber < 8 + structure.atoms_number; ++currLineNumber) {
    if (infile.ignore(std::numeric_limits<std::streamsize>::max(),
                      infile.widen('\n'))) {
      // just skipping the line. POSCAR must be read and checked seperately
    } else {
      log->error(file_name + " could not be read properly");
      return {};
    }
  }
  log->trace("Started reading " + file_name);
  std::getline(infile, temp_line);
  arma::urowvec3 grid;
  infile >> grid;
  arma::cube rawdata_cube(as_size(grid));
  infile >> rawdata_cube;
  std::getline(infile, temp_line);

  return rawdata_cube;
}

void supercell::write_CHGPOT(const std::string &type,
                             const std::string &file_name) const {
  auto log = spdlog::get("loggers");
  log->trace("Started writing " + file_name);

  write_POSCAR(file_name);
  std::ofstream out_file;
  out_file.open(file_name, std::ofstream::app);

  arma::cube CHGPOT;
  if (type == "CHGCAR") {
    CHGPOT = charge;
  } else if (type == "LOCPOT") {
    CHGPOT = potential;
  }
  out_file << '\n' << SizeVec(CHGPOT) << '\n';
  out_file << std::fixed << std::showpos << std::scientific
           << std::setprecision(10);
  for (arma::uword i = 0; i < CHGPOT.n_elem; ++i) {
    out_file << CHGPOT(i) << " ";
    if (i % 5 == 4) {
      out_file << '\n';
    }
  }
  out_file.close();
}

void supercell::write_CHGCAR(const std::string &file_name) const {
  write_CHGPOT("CHGCAR", file_name);
}

void supercell::write_LOCPOT(const std::string &file_name) const {
  write_CHGPOT("LOCPOT", file_name);
}
void write_planar_avg(const arma::cube &potential_data,
                      const arma::cube &charge_data, const std::string &id,
                      const arma::rowvec3 &coordinate_vectors,
                      const int direction) {
  auto log = spdlog::get("loggers");
  unsigned int direction_first = 0;
  unsigned int direction_last = 2;

  if (direction != -1) {
    direction_first = direction;
    direction_last = direction;
  }

  for (unsigned int dir = direction_first; dir <= direction_last; ++dir) {
    auto avg_pot = planar_average(dir, potential_data);
    const auto avg_chg = planar_average(dir, charge_data);
    const auto pot_normalization =
        static_cast<double>(avg_pot.n_elem) / potential_data.n_elem;
    avg_pot *= pot_normalization;

    arma::vec pot_coordinates = arma::linspace<arma::vec>(
        0, coordinate_vectors(dir) / ang_to_bohr, avg_pot.n_elem + 1);
    arma::vec chg_coordinates = arma::linspace<arma::vec>(
        0, coordinate_vectors(dir) / ang_to_bohr, avg_chg.n_elem + 1);
    arma::mat pot_data = arma::zeros(avg_pot.n_elem, 2);
    for (arma::uword i = 0; i < avg_pot.n_elem; ++i) {
      pot_data.row(i) = arma::rowvec{pot_coordinates(i), avg_pot(i)};
    }
    arma::mat chg_data = arma::zeros(avg_chg.n_elem, 2);
    for (arma::uword i = 0; i < avg_pot.n_elem; ++i) {
      chg_data.row(i) = arma::rowvec{chg_coordinates(i), avg_chg(i)};
    }
    const std::string pot_filename = "slabcc_" + id + int2xyz(dir) + "POT.dat";
    const std::string chg_filename = "slabcc_" + id + int2xyz(dir) + "CHG.dat";
    log->trace("Started writing {}", pot_filename);
    pot_data.save(pot_filename, arma::raw_ascii);
    log->trace("Started writing {}", chg_filename);
    chg_data.save(chg_filename, arma::raw_ascii);
  }
}

void check_slabcc_compatiblity(const supercell &Neutral_supercell,
                               const supercell &Charged_supercell) {

  auto log = spdlog::get("loggers");

  // equal size
  if (!approx_equal(Neutral_supercell.cell_vectors * Neutral_supercell.scaling,
                    Charged_supercell.cell_vectors * Charged_supercell.scaling,
                    "reldiff", 0.001)) {
    log->debug(
        "Neutral supercell vectors: " +
        to_string(Neutral_supercell.cell_vectors * Neutral_supercell.scaling));
    log->debug(
        "Charged supercell vectors: " +
        to_string(Charged_supercell.cell_vectors * Charged_supercell.scaling));
    log->critical("Cell vectors of the input files does not match!");
    finalize_loggers();
    exit(1);
  }

  // orthogonal
  auto cell = Neutral_supercell.cell_vectors;
  cell.each_col([](arma::vec &vec) { vec /= arma::norm(vec); });
  if (!arma::approx_equal(cell.t() * cell, arma::eye(3, 3), "reldiff", 0.001)) {
    log->debug("Cell vectors: " + to_string(Neutral_supercell.cell_vectors));
    log->debug("Cell vectors basis: " + to_string(cell));
    log->debug("Orthogonality criteria: " + to_string(cell.t() * cell));
    log->critical("Supercell vectors are not orthogonal!");
    finalize_loggers();
    exit(1);
  }

  // equal grid
  const arma::SizeCube input_grid = arma::size(Neutral_supercell.potential);
  if ((input_grid != arma::size(Charged_supercell.potential)) ||
      (input_grid != arma::size(Charged_supercell.charge)) ||
      (input_grid != arma::size(Neutral_supercell.charge))) {
    log->debug("Neutral CHGCAR grid: " +
               to_string(arma::size(Neutral_supercell.charge)));
    log->debug("Neutral LOCPOT grid: " +
               to_string(arma::size(Neutral_supercell.potential)));
    log->debug("Charged CHGCAR grid: " +
               to_string(arma::size(Charged_supercell.charge)));
    log->debug("Charged LOCPOT grid: " +
               to_string(arma::size(Charged_supercell.potential)));
    log->critical(
        "Grid size of the data in CHGCAR/LOCPOT files does not match!");
    finalize_loggers();
    exit(1);
  }

  log->trace("All files are loaded and cross-checked!");
}
