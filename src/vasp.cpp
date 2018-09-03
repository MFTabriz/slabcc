// Copyright (c) 2018, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "vasp.hpp"

void write_POSCAR(const supercell& structure, const string& file_name) {
	ofstream out_file;
	out_file.open(file_name);
	out_file << structure.label << endl;
	out_file << "   " << fixed << setprecision(15) << structure.scaling << endl;
	structure.cell_vectors.t().raw_print(out_file);
	vector<int> defs = structure.atoms.definition_order;
	auto end = unique(defs.begin(), defs.end(), [](const int& l, const int& r) noexcept { return  l == r; });
	defs.erase(end, defs.end());
	for (const int& i : defs) {
		auto it = find(structure.atoms.definition_order.begin(), structure.atoms.definition_order.end(), i);
		const auto pos = distance(structure.atoms.definition_order.begin(), it);
		out_file << " " << structure.atoms.type.at(pos);
	}
	out_file << endl;

	int counter = 0;
	int counted = 0;
	for (int k = 0; k < defs.size(); ++k) {
		while (structure.atoms.definition_order.at(counter) == defs.at(k)) {
			++counter;
			if (counter == structure.atoms_number) break;
		}
		out_file << " " << counter - counted;
		counted = counter;
	}
	out_file << endl;

	if (structure.selective_dynamics) {
		out_file << "Selective dynamics" << endl;
	}
	out_file << structure.coordination_system << endl;

	for (uword number = 0; number < structure.atoms_number; ++number) {
		out_file << arma::rowvec(structure.atoms.position.row(number));
		if (structure.selective_dynamics) {
			for (int qq2 = 0; qq2 < 3; ++qq2) {
				if (structure.atoms.constrains.at(number).at(qq2)) {
					out_file << " F";
				}
				else {
					out_file << " T";
				}
			}
		}
		out_file << endl;
	}

	out_file.close();
}


rowvec3 direct_cord(const supercell& structure, const rowvec3& cartesians) {
	mat33 cell_vectors = structure.cell_vectors / structure.scaling;
	rowvec3 direct_coords = solve(cell_vectors.t(), cartesians.t()).t();
	return direct_coords;
}

void normalize_positions(supercell& structure) {
	if (structure.coordination_system == "direct") {
		structure.atoms.position = fmod_p(structure.atoms.position, 1);
	}
	else {
		cout << "WARNING: positions in cartesian coordinates cannot be normalized!" << endl;
	}
}

supercell read_POSCAR(const string& file_name) {
	supercell structure;
	ifstream infile;
	string temp_line;
	infile.open(file_name);
	if (!infile) {
		cerr << "While opening the POSCAR file an error is encountered!" << endl;
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
	for (auto i = 0; i < temp_types.size(); ++i) {
		for (int a = 0; a < temp_numbers.at(i); ++a) {
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
		cerr << "ERROR: Unknown coordination system type!" << endl;
		return structure;
	}

	for (uword i = 0; i < structure.atoms_number; ++i) {
		string temp_constrain;
		infile >> structure.atoms.position.row(i);
		if (structure.selective_dynamics) {
			for (int qq2 = 0; qq2 < 3; ++qq2) {
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
	ifstream infile;
	string temp_line;
	infile.open(file_name);
	if (!infile) {
		cerr << "ERROR: file not found:" << file_name << endl;
		return {};
	}
	supercell structure = read_POSCAR(file_name);
	for (auto currLineNumber = 0; currLineNumber < 8 + structure.atoms_number; ++currLineNumber) {
		if (infile.ignore(numeric_limits<streamsize>::max(), infile.widen('\n'))) {
			//just skipping the line. POSCAR must be read and checked seperately
		}
		else {
			cerr << "ERROR: Could not read file properly: " << file_name << endl;
			return {};
		}
	}
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
	supercell structure_info = structure;

	write_POSCAR(structure_info, file_name);
	ofstream out_file;
	out_file.open(file_name, ofstream::app);

	cube CHGPOT;
	if (type == "CHGCAR") {
		CHGPOT = structure.charge;
	}
	else if (type == "LOCPOT") {
		CHGPOT = structure.potential;
	}
	out_file << endl << SizeVec(CHGPOT) << endl;
	out_file << fixed << showpos << scientific << setprecision(10);
	for (uword i = 0; i < CHGPOT.n_elem; ++i) {
		out_file << CHGPOT(i) << " ";
		if (i % 5 == 4) {
			out_file << endl;
		}
	}
	out_file.close();
}

void write_planar_avg(const cube& potential_data, const cube& charge_data, const string& id) {

	vector<double> AV1 = planar_average(0, potential_data);
	vector<double> AV2 = planar_average(1, potential_data);
	vector<double> AV3 = planar_average(2, potential_data);
	for (auto &elem : AV1) {
		elem /= potential_data.n_cols * potential_data.n_slices;
	}
	for (auto &elem : AV2) {
		elem /= potential_data.n_rows * potential_data.n_slices;
	}
	for (auto &elem : AV3) {
		elem /= potential_data.n_rows * potential_data.n_cols;
	}


	write_vec2file(AV1, "slabcc_" + id + "XPOT.dat");
	write_vec2file(AV2, "slabcc_" + id + "YPOT.dat");
	write_vec2file(AV3, "slabcc_" + id + "ZPOT.dat");

	AV1 = planar_average(0, charge_data);
	AV2 = planar_average(1, charge_data);
	AV3 = planar_average(2, charge_data);

	write_vec2file(AV1, "slabcc_" + id + "XCHG.dat");
	write_vec2file(AV2, "slabcc_" + id + "YCHG.dat");
	write_vec2file(AV3, "slabcc_" + id + "ZCHG.dat");
}
