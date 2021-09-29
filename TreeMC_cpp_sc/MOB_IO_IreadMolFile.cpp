#pragma once
#include "MOB_IO_IreadMolFile.h"

namespace MOB_IO_IreadMolFile {

	MolData parse_mol_file(std::string mol_file_name, bool include_H) {
		MolData mol_data;

		// open mol file and check it exists
		std::ifstream mol_file;
		mol_file.open(mol_file_name);
		if (!mol_file.good()) {
			mol_data.good = false;
			return mol_data;
		}

		// initialise stuff required to iterate through file and fill atoms/bonds
		std::string current_line;
		auto line{ 1 }, current_atom{ 0 }, current_bond{ 0 }, atom_block_end{ 0 }, bond_block_end{ 0 }, num_atoms{ 0 }, num_bonds{ 0 };
		std::map<int, std::string> atoms;
		std::map<int, std::vector<int>> bonds;
		std::vector<std::string> atom_vec;
		std::vector<std::vector<int>> bond_vec;

		while (std::getline(mol_file, current_line)) {

			// line 4 contains the number of atoms and bonds. Parse this, and work out the end position of the atom and bond blocks
			if (line == 4) {
				auto num_atom_string = current_line.substr(0, 3);
				auto num_bond_string = current_line.substr(3, 3);
				num_atoms = atoi(num_atom_string.c_str());
				num_bonds = atoi(num_bond_string.c_str());
				atom_block_end = line + num_atoms;
				bond_block_end = atom_block_end + num_bonds;
				/*std::cout << "num atoms " << num_atoms << ", num bonds " << num_bonds << std::endl;*/

			}
			else if (line > 4 && line <= atom_block_end) {
				// iterate through the atom block

				std::string current_atom_type;

				// atom type may be 1 or 2 characters. If 2, the second postion is a space
				if (current_line.substr(32, 1) == " ") {
					current_atom_type = current_line.substr(31, 1);
				}
				else {
					current_atom_type = current_line.substr(31, 2);
				}
				atoms[current_atom] = current_atom_type;
				current_atom++;
			}

			else if ((line > atom_block_end) && (line <= bond_block_end)) // go through the bond block
			{
				// atoms at position 1 and 4, bond type at position 8
				int bond_type = atoi(current_line.substr(8, 1).c_str());
				int atom1 = std::stoi(current_line.substr(0, 3)) - 1; // -1 due to zero-indexing
				int atom2 = std::stoi(current_line.substr(3, 3)) - 1; // -1 due to zero-indexing
				bonds[current_bond] = { atom1, atom2, bond_type };

				current_bond++;
			}

			line += 1;
		}


		if (include_H) {
			// no more work to do if including H atoms, convert to vectors and return
			// note - if there are no H atoms in the mol file in the first place, they won't be added
			map_to_vector(atoms, atom_vec);
			map_to_vector(bonds, bond_vec);
			return MolData(atom_vec, bond_vec);

		}
		else {
			// make new maps without the H and with the indices renumbered
			std::map<int, std::string> renumbered_atoms;
			std::map<int, std::vector<int>>renumbered_bonds;
			std::map<int, int> index_map; // this is so the bond indices can be correctly renumbered
			int new_atom_index{ 0 }, new_bond_index{ 0 };

			for (auto& item : atoms) {
				// only keep non-H atoms
				if (item.second != "H") {
					renumbered_atoms[new_atom_index] = item.second;
					index_map[item.first] = new_atom_index;
					new_atom_index++;
				}
			}
			for (auto& item : bonds) {
				// only keep bonds that don't contain H
				auto a_1_index = item.second[0];
				auto a_2_index = item.second[1];
				auto b_type = item.second[2];
				auto a_1 = atoms[a_1_index];
				auto a_2 = atoms[a_2_index];
				if (a_1 != "H" && a_2 != "H") {
					renumbered_bonds[new_bond_index] = { index_map[a_1_index], index_map[a_2_index], b_type };
					new_bond_index++;
				}
			}

			// convert the maps to the required vectors and return
			map_to_vector(renumbered_atoms, atom_vec);
			map_to_vector(renumbered_bonds, bond_vec);
			return MolData(atom_vec, bond_vec);


		}
	}



	//void test_print(std::map<int, std::string> atoms, std::map<int, std::vector<int>> bonds) {
	//	std::cout << "Atoms" << std::endl;
	//	for (auto& item : atoms) {
	//		std::cout << item.first << ": " << item.second << std::endl;
	//	}
	//	std::cout << "Bonds" << std::endl;
	//	for (auto& item : bonds) {
	//		std::cout << item.first << ", " << "[";
	//		for (auto& b : item.second) {
	//			std::cout << b << ", ";
	//		}
	//		std::cout << "]" << std::endl;
	//	}
	//}

	//void test_print_mol(MOB_IO_IreadMolFile::MolData m) {
	//	std::cout << "atoms" << std::endl;
	//	for (auto& atom : m.atoms) {
	//		std::cout << atom << std::endl;
	//	}std::cout << "bonds" << std::endl;
	//	for (auto& bond : m.bonds) {
	//		for (auto& b : bond) {
	//			std::cout << b << " ";
	//		}
	//		std::cout << std::endl;
	//	}
	//}

	template <typename M, typename V>
	void map_to_vector(const M& m, V& v) {
		for (auto& item : m) {
			v.push_back(item.second);
		}

	}
}
