// -- PCTree associated function --
#include "MOB.h"
#include <iostream>


/* =========================================== */
int MOL_BOND::addMol(strType molFile) {
	MOL_BOND obj{};
	std::vector<MOL_BOND> objList{};
	int intMsg = obj.importData(molFile, objList);
	if (intMsg != 1) {
		throw " ! Error: Mol-file contains 0 or more than 1 molecule.";
	}
	else {
		obj = objList[0];
	}

	// find the max ID +1 (set to be the starting ID next)
	int atomIDmax1 = -1;
	for (auto& atom : Atoms)
		if (atomIDmax1 < atom.first)
			atomIDmax1 = atom.first;
	atomIDmax1++;
	int bondIDmax1 = -1;
	for (auto& bond : Bonds)
		if (bondIDmax1 < bond.first)
			bondIDmax1 = bond.first;
	bondIDmax1++;

	// make the new IDs for Atoms and Bonds
	for (auto bond = obj.Bonds.begin(); bond != obj.Bonds.end(); bond++) {
		bond->second.Atom1id += atomIDmax1;
		bond->second.Atom2id += atomIDmax1;
		Bonds.insert({ bondIDmax1 + bond->first, bond->second });
	}
	for (auto atom = obj.Atoms.begin(); atom != obj.Atoms.end(); atom++) {
		for (int& i : atom->second.BondsLinked) {
			i += bondIDmax1;
		}
		for (auto& vec2d : atom->second.toCut) {
			for (intVec& vec1d : vec2d) {
				for (int& i : vec1d) {
					i += bondIDmax1;
				}
			}
		}
		Atoms.insert({ atomIDmax1 + atom->first, atom->second });
	}


	// change AtomsHub;
	for (int& i : obj.AtomsHub) {
		AtomsHub.push_back(i + atomIDmax1);
	}

	// change the other to undefined
	inchi = "";
	nAtomsIncludeH = -1;

	return bondIDmax1;
}


/* =========================================== */
std::vector<PATHWAY_ASSEMBLY> MOL_BOND::collectPaths(
	FRAGHIST& fraghist, int Nmol, 
	int NtryPath, int Ntry_freeze) {

	PATHWAY_ASSEMBLY PAempty{};
	PAempty.trivial_index = static_cast<int>(Bonds.size()) - Nmol;
	std::vector<PATHWAY_ASSEMBLY> PAvec{};
	int indexMin = INT_MAX;

	for (int i = 0; i < NtryPath; ++i) {
		
		PATHSINGLE path = get1pathway(fraghist, Nmol, Ntry_freeze);
		if (indexMin > path.PAindex) {
			indexMin = path.PAindex;
		}
		std::cout << i << ":" << path.PAindex << "(" << indexMin << "), ";

		if (PAvec.size() <= path.PAindex) {
			for (size_t j = PAvec.size(); j <= path.PAindex; ++j) {
				PAvec.push_back(PAempty);
			}
		}
		PATHWAY_ASSEMBLY& thePA = PAvec[path.PAindex];
		if (thePA.paths.empty()) {
			thePA.paths.push_back(path);
		}
		else {
			bool isPathNew = true;
			for (size_t j = 0; j < thePA.paths.size(); ++j) {
				if (isPathEqual(path, thePA.paths[j])) {
					isPathNew = false;
					break;
				}
			}
			if (isPathNew) { // if path is new, then add it to PA.pahts
				thePA.paths.push_back(path);
			}
		}
	}
	std::cout << std::endl;
	return PAvec;
}
