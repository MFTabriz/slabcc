// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "vasp.hpp"

void write_POSCAR(const supercell& structure, const string& file_name) {
	ofstream out_file;
	out_file.open(file_name);
	out_file << structure.label << '\n';
	out_file << "   " << fixed << showpos << setprecision(16) << structure.scaling << '\n';
	structure.cell_vectors.t().raw_print(out_file);
	vector<int> defs = structure.atoms.definition_order;
	const auto end = unique(defs.begin(), defs.end(), [](const int& l, const int& r) noexcept { return  l == r; });
	defs.erase(end, defs.end());
	for (const int& i : defs) {
		const auto it = find(structure.atoms.definition_order.begin(), structure.atoms.definition_order.end(), i);
		const auto pos = distance(structure.atoms.definition_order.begin(), it);
		out_file << " " << structure.atoms.type.at(pos);
	}
	out_file << '\n' << noshowpos;

	uword counter = 0;
	uword counted = 0;
	for (size_t k = 0; k < defs.size(); ++k) {
		while (structure.atoms.definition_order.at(counter) == defs.at(k)) {
			++counter;
			if (counter == structure.atoms_number) break;
		}
		out_file << " " << counter - counted;
		counted = counter;
	}
	out_file << '\n';

	if (structure.selective_dynamics) {
		out_file << "Selective dynamics\n";
	}
	out_file << structure.coordination_system << '\n';

	for (uword number = 0; number < structure.atoms_number; ++number) {
		out_file << " " << structure.atoms.position.row(number);
		if (structure.selective_dynamics) {
			for (auto qq2 = 0; qq2 < 3; ++qq2) {
				if (structure.atoms.constrains.at(number).at(qq2)) {
					out_file << " F";
				}
				else {
					out_file << " T";
				}
			}
		}
		out_file << '\n';
	}

	out_file.close();
}

rowvec3 direct_cord(const supercell& structure, const rowvec3& cartesians) {
	const mat33 cell_vectors = structure.cell_vectors / structure.scaling;
	const rowvec3 direct_coords = solve(cell_vectors.t(), cartesians.t()).t();
	return direct_coords;
}

void normalize_positions(supercell& structure) {
	auto log = spdlog::get("loggers");
	if (structure.coordination_system == "direct") {
		structure.atoms.position = fmod_p(structure.atoms.position, 1);
	}
	else {
		//I'm lying here! > It can be converted to direct coordinates, normalized and converted back.
		log->warn("positions in the Cartesian coordinates cannot be normalized!");
	}
}

supercell read_POSCAR(const string& file_name) {
	auto log = spdlog::get("loggers");
	supercell structure;
	ifstream infile;
	string temp_line;
	infile.open(file_name);
	if (!infile) {
		log->critical("Could not open the "+ file_name);
		return structure;
	}
	getline(infile, structure.label);
	getline(infile, temp_line);
	structure.scaling = stod(temp_line);
	infile >> structure.cell_vectors;

	getline(infile, temp_line);
	getline(infile, temp_line);

	string buf;
	stringstream ss(temp_line);
	vector <string> temp_types;
	vector <int> temp_numbers;

	while (ss >> buf) {
		temp_types.push_back(buf);
	}
	getline(infile, temp_line);
	stringstream ss2(temp_line);
	while (ss2 >> buf) {
		temp_numbers.push_back(stoi(buf));
	}

	int atom_counter = 0;
	for (size_t i = 0; i < temp_types.size(); ++i) {
		for (auto a = 0; a < temp_numbers.at(i); ++a) {
			structure.atoms.type.push_back(temp_types.at(i));
			structure.atoms.definition_order.push_back(i);
			++atom_counter;
		}
	}
	structure.atoms_number = atom_counter;
	structure.atoms.position.set_size(atom_counter, 3);
	structure.atoms.constrains.resize(atom_counter, vector<bool>(3, 0));
	getline(infile, temp_line);
	if (tolower(temp_line.at(0)) == 's') {
		structure.selective_dynamics = true;
		getline(infile, temp_line);
	}
	else {
		structure.selective_dynamics = false;
	}

	if ((tolower(temp_line.at(0)) == 'k') || (tolower(temp_line.at(0)) == 'c')) {
		structure.coordination_system = "cartesian";
	}
	else if (tolower(temp_line.at(0)) == 'd') {
		structure.coordination_system = "direct"; //in VASP it means relative!!
	}
	else {
		log->critical("Unknown coordination system type in "+ file_name);
		return structure;
	}

	for (uword i = 0; i < structure.atoms_number; ++i) {
		string temp_constrain;
		infile >> structure.atoms.position.row(i);
		if (structure.selective_dynamics) {
			for (auto qq2 = 0; qq2 < 3; ++qq2) {
				infile >> temp_constrain;
				structure.atoms.constrains.at(i).at(qq2) = (tolower(temp_constrain.at(0)) == 't') ? false : true;
			}
		}
		else {
			infile.ignore(numeric_limits<streamsize>::max(), '\n');
		}
	}

	normalize_positions(structure);
	return structure;

}

cube read_CHGPOT(const string& file_name) {
	auto log = spdlog::get("loggers");
	ifstream infile;
	string temp_line;
	infile.open(file_name);
	if (!infile) {
		log->critical("File not found: " + file_name);
		return {};
	}
	const supercell structure = read_POSCAR(file_name);
	for (uword currLineNumber = 0; currLineNumber < 8 + structure.atoms_number; ++currLineNumber) {
		if (infile.ignore(numeric_limits<streamsize>::max(), infile.widen('\n'))) {
			//just skipping the line. POSCAR must be read and checked seperately
		}
		else {
			log->error(file_name+ " could not be read properly");
			return {};
		}
	}
	log->trace("Started reading "+ file_name);
	getline(infile, temp_line);
	urowvec3 grid;
	infile >> grid;
	cube rawdata_cube(as_size(grid));
	infile >> rawdata_cube;
	getline(infile, temp_line);

	return rawdata_cube;
}

void shift_structure(supercell& structure, const rowvec3& relative_shift) {

	rowvec3 pos_shift = relative_shift;

	// if our atomic positions are in cartesian coordinates,
	// find the equivalent cartesian vector
	if (structure.coordination_system == "cartesian") {
		pos_shift = direct_cord(structure, pos_shift);
	}

	structure.atoms.position.each_row() += pos_shift;

	if (structure.coordination_system == "direct") {
		normalize_positions(structure);
	}

	structure.charge = shift(structure.charge, relative_shift);
	structure.potential = shift(structure.potential, relative_shift);

}

void write_CHGPOT(const string& type, const string& file_name, const supercell& structure) {
	auto log = spdlog::get("loggers");
	log->trace("Started writing " + file_name);

	write_POSCAR(structure, file_name);
	ofstream out_file;
	out_file.open(file_name, ofstream::app);

	cube CHGPOT;
	if (type == "CHGCAR") {
		CHGPOT = structure.charge;
	}
	else if (type == "LOCPOT") {
		CHGPOT = structure.potential;
	}
	out_file << '\n' << SizeVec(CHGPOT) << '\n';
	out_file << fixed << showpos << scientific << setprecision(10);
	for (uword i = 0; i < CHGPOT.n_elem; ++i) {
		out_file << CHGPOT(i) << " ";
		if (i % 5 == 4) {
			out_file << '\n';
		}
	}
	out_file.close();
}

void write_planar_avg(const cube& potential_data, const cube& charge_data, const string& id, const int direction) {
	unsigned int direction_first = 0;
	unsigned int direction_last = 2;

	if (direction != -1) {
		direction_first = direction;
		direction_last = direction;
	}

	for (unsigned int dir = direction_first; dir <= direction_last; ++dir) {
		vector<double> avg_pot = planar_average(dir, potential_data);
		const vector<double> avg_chg = planar_average(dir, charge_data);
		const auto pot_normalization = static_cast<double>(avg_pot.size()) / potential_data.n_elem;
		for (auto &elem : avg_pot) {
			elem *= pot_normalization;
		}
		write_vec2file(avg_pot, "slabcc_" + id + int2xyz(dir) + "POT.dat");
		write_vec2file(avg_chg, "slabcc_" + id + int2xyz(dir) + "CHG.dat");
	}
}



void check_slabcc_compatiblity(const supercell& Neutral_supercell, const supercell& Charged_supercell) {

	auto log = spdlog::get("loggers");

	// equal size
	if (!approx_equal(Neutral_supercell.cell_vectors * Neutral_supercell.scaling,
		Charged_supercell.cell_vectors * Charged_supercell.scaling, "reldiff", 0.001)) {
		log->debug("Neutral supercell vectors: " + to_string(Neutral_supercell.cell_vectors * Neutral_supercell.scaling));
		log->debug("Charged supercell vectors: " + to_string(Charged_supercell.cell_vectors * Charged_supercell.scaling));
		log->critical("Cell vectors of the input files does not match!");
		finalize_loggers();
		exit(1);
	}

	//orthogonal
	auto cell = Neutral_supercell.cell_vectors;
	cell.each_col([](vec & vec) { vec /= norm(vec); });
	if (!approx_equal(cell.t() * cell, eye(3, 3), "reldiff", 0.001)) {
		log->debug("Cell vectors: " + to_string(Neutral_supercell.cell_vectors));
		log->debug("Cell vectors basis: " + to_string(cell));
		log->debug("Orthogonality criteria: " + to_string(cell.t() * cell));
		log->critical("Supercell vectors are not orthogonal!");
		finalize_loggers();
		exit(1);
	}

	// equal grid
	const SizeCube input_grid = arma::size(Neutral_supercell.potential);
	if ((input_grid != arma::size(Charged_supercell.potential)) ||
		(input_grid != arma::size(Charged_supercell.charge)) ||
		(input_grid != arma::size(Neutral_supercell.charge))) {
		log->debug("Neutral CHGCAR grid: " + to_string(arma::size(Neutral_supercell.charge)));
		log->debug("Neutral LOCPOT grid: " + to_string(arma::size(Neutral_supercell.potential)));
		log->debug("Charged CHGCAR grid: " + to_string(arma::size(Charged_supercell.charge)));
		log->debug("Charged LOCPOT grid: " + to_string(arma::size(Charged_supercell.potential)));
		log->critical("Data grids of the CHGCAR or LOCPOT files does not match!");
		finalize_loggers();
		exit(1);
	}

	log->trace("All files are loaded and cross-checked!");
}
