/* IO_I_readMolFile is a namespace which contains classes and functions to read mol-file.txt.
   It reads mol-file.txt and returns a MolData class,
   which can be used by function getMOL_from() [in func1_InOut.cpp]
   to transform it into a MOL instance.
*/
#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>


namespace MOB_IO_IreadMolFile {

	class MolData {
	public:
		MolData() {};
		MolData(std::vector<std::string> a, std::vector<std::vector<int>> b, bool g = true) : atoms(a), bonds(b), good(g) {}

		std::vector<std::string> atoms;
		std::vector<std::vector<int>> bonds;
		bool good = false;
	};



	MolData parse_mol_file(std::string mol_file_name, bool include_H = false);
	//void test_print(std::map<int, std::string> atoms, std::map<int, std::vector<int>> bonds);

	template <typename M, typename V>
	void map_to_vector(const M& m, V& v);
	//void test_print_mol(MolData m);
}