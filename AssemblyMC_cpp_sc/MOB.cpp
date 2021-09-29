#include "MOB.h"
#include "MOB_IO_IreadMolFile.h"
#include "MOB_IO_IconvertToMOB.h"
#include "MOB_Inchi.h"
#include <stack>
#include <iostream>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <sstream>


/* ============== associated functions ================= */
// find_x_in_V() is intermediate function that checks whether 
// integer x is in the vector V
bool find_x_in_V(int x, intVec& V) {
	for (size_t i = 0; i < V.size(); ++i) {
		if (V[i] == x) {
			return true;
		}
	}
	return false;
}

/* =========================================== */
// correspond to getToCut()
size_t getToCutInt(size_t n) {
	switch (n) {
	case 1:
		return 1;
	case 2:
		return 2;
	case 3:
		return 5;
	case 4:
		return 15;
	case 5:
		return 52;
	default:
		return 0;
	}
}

/* =========================================== */
strType getIdentifierStr(intVec& identifierVec) {
	strType identifier{};
	std::sort(identifierVec.begin(), identifierVec.end());
	for (int& i : identifierVec) {
		identifier += (std::to_string(i) + ',');
	}
	return identifier;
}
intVec retrieveBondsID(strType identifierStr) {
	intVec bondsID{};
	strType token;
	std::istringstream tokenStream(identifierStr);
	while (std::getline(tokenStream, token, ',')) {
		bondsID.push_back(std::stoi(token));
	}
	return bondsID;
}

/* =========================================== */
// check if vec is inside vec2D
// note it's better that item in vec2D has the same length as vec
bool isVecInVec2D(intVec& vec, intVec2D& vec2D) {
	for (intVec& thisVec : vec2D) {
		if (thisVec == vec) {
			return true;
		}
	}
	return false;
}
// check if child vector is included in mother vector
bool isEmbracing(intVec& mother, intVec& child) {
	size_t j = 0;
	for (int& i : child) {
		bool found = false;
		while (j < mother.size()) {
			if (i == mother[j]) {
				found = true;
				break;
			}
			++j;
		}
		if (!found)
			return false;
		++j;
	}
	return true;
}

/* =========================================== */
// find fragID's corresponding index in x0 that stores the list of inchis
size_t findIDsIndexInX0(intVec& fragID,
	std::vector<strType>& x0,
	std::unordered_map<strType, intVec2D>& xID0) {
	for (size_t i = 0; i < x0.size(); ++i) {
		for (intVec& id : xID0[x0[i]]) {
			if (id == fragID) {
				return i;
			}
		}
	}
	return 0;
}

/* =========================================== */
// check whether two paths are equivalent
bool isPathEqual(PATHSINGLE& path1, PATHSINGLE& path2) {
	if (path1.PAindex == path2.PAindex
		&& path1.dupN.size() == path2.dupN.size()) {
		for (size_t i = 0; i < path1.dup.size(); ++i) {
			bool found = false;
			for (size_t j = 0; j < path2.dup.size(); ++j) {
				if (path1.dupN[i] == path2.dupN[j] &&
					path1.dup[i] == path2.dup[j]) {
					found = true;
					break;
				}
			}
			if (!found) {
				return false;
			}
		}
		return true;
	}
	else {
		return false;
	}
}

/* =========================================== */
// Get the total number of atoms including H atoms from InChI string.
int getTotalNumAtoms(MOL_BOND& mol, strType& inc) {
	bool foundH = false;
	auto it = inc.begin();
	while (it != inc.end()) {
		if (*it == 'H') {
			foundH = true;
			break;
		}
		it++;
	}
	// After this, the letter 'H' is found (foundH is true now), 
	// and "it" points to the letter after "H".

	if (foundH) {
		it++;
		strType temp = "";
		// read the integer after "H" until a letter,
		// so the number of H atoms is get in "temp" as a string.
		while (it != inc.end()) {
			if (*it == '0' || *it == '1' || *it == '2' || *it == '3' ||
				*it == '4' || *it == '5' || *it == '6' || *it == '7' ||
				*it == '8' || *it == '9') {
				temp.push_back(*it);
				it++;
			}
			else { break; }
		}
		if (temp == "") {
			temp = '0';
		}
		int nAtoms = static_cast<int>(mol.Atoms.size() + stoi(temp));
		return nAtoms;
	}
	else {
		// The letter "H" is not found in the InChI string.
		return static_cast<int>(mol.Atoms.size());
	}
}

/* =========================================== */
// to check tryThis fragment has been filled or not,
// or links several fragments
// the first 2 are input parameters, the last variable can be changed
int tryGunshot(intVec& tryThis, intVec2D& target, intVec& linkedTo) {
	int iFound = -1, jFound = -1, kFound = -1;
	bool Found = false;
	for (int k = 0; k < tryThis.size(); ++k) {
		for (int i = 0; i < target.size(); ++i) {
			for (int j = 0; j < target[i].size(); ++j) {
				if (tryThis[k] == target[i][j]) {
					iFound = i;
					jFound = j;
					kFound = k;
					linkedTo[k] = i;
					Found = true;
					break;
				}
			}
			if (Found)
				break;
		}
		if (Found)
			break;
	}
	if (!Found)
		return 1; // tryThis is available

	// check if tryThis is totally covered (treat as a special case)
	if (kFound == 0) {
		bool case2 = true;
		for (size_t k = 1; k < tryThis.size(); ++k) {
			Found = false;
			for (int j = jFound + 1; j < target[iFound].size(); ++j) {
				if (tryThis[k] == target[iFound][j]) {
					Found = true;
					break;
				}
			}
			if (!Found) {
				case2 = false;
				break;
			}
		}
		if (case2)
			return 2; // tryThis is totally included in an already-filled fragment
	}

	// other case is case 3
	for (int k = kFound + 1; k < tryThis.size(); ++k) {
		for (int i = 0; i < target.size(); ++i) {
			for (size_t j = 0; j < target[i].size(); ++j) {
				if (tryThis[k] == target[i][j]) {
					linkedTo[k] = i;
				}
			}
		}
	}
	return 3;
}

/* =========================================== */
std::vector<intVec2D> getToCut(intVec& BondsLinked) {
	int a, b, c, d, e;
	switch (BondsLinked.size()) {
	case 2:
		a = BondsLinked[0];
		b = BondsLinked[1];
		return std::vector<intVec2D>{ // = 2 ways
			intVec2D{ intVec{ a,b } },
				intVec2D{ intVec{ a }, intVec{ b } },
		};
	case 3:
		a = BondsLinked[0];
		b = BondsLinked[1];
		c = BondsLinked[2];
		return std::vector<intVec2D>{ // = 5 ways
			intVec2D{ intVec{ a,b,c } },
				intVec2D{ intVec{ a }, intVec{ b,c } },
				intVec2D{ intVec{ b }, intVec{ a,c } },
				intVec2D{ intVec{ c }, intVec{ a,b } },
				intVec2D{ intVec{ a }, intVec{ b }, intVec{ c } }
		};
	case 4:
		a = BondsLinked[0];
		b = BondsLinked[1];
		c = BondsLinked[2];
		d = BondsLinked[3];
		return std::vector<intVec2D>{ // = 15 ways
			intVec2D{ intVec{ a,b,c,d } }, // 4
				intVec2D{ intVec{ a }, intVec{ b,c,d } },// 1+3
				intVec2D{ intVec{ b }, intVec{ a,c,d } },
				intVec2D{ intVec{ c }, intVec{ a,b,d } },
				intVec2D{ intVec{ d }, intVec{ a,b,c } },
				intVec2D{ intVec{ a,b }, intVec{ c,d } },// 2+2
				intVec2D{ intVec{ a,c }, intVec{ b,d } },
				intVec2D{ intVec{ a,d }, intVec{ b,c } },
				intVec2D{ intVec{ c }, intVec{ d }, intVec{ a,b } },// 1+1+2
				intVec2D{ intVec{ b }, intVec{ d }, intVec{ a,c } },
				intVec2D{ intVec{ b }, intVec{ c }, intVec{ a,d } },
				intVec2D{ intVec{ a }, intVec{ d }, intVec{ b,c } },
				intVec2D{ intVec{ a }, intVec{ c }, intVec{ b,d } },
				intVec2D{ intVec{ a }, intVec{ b }, intVec{ c,d } },
				intVec2D{ intVec{ a }, intVec{ b }, intVec{ c }, intVec{ d } }// 1+1+1+1
		};
	case 5:
		a = BondsLinked[0];
		b = BondsLinked[1];
		c = BondsLinked[2];
		d = BondsLinked[3];
		e = BondsLinked[4];
		return std::vector<intVec2D>{ // = 52 ways
			intVec2D{ intVec{ a,b,c,d,e } }, // 5
				intVec2D{ intVec{ a }, intVec{ b,c,d,e } },// 1+4
				intVec2D{ intVec{ b }, intVec{ a,c,d,e } },
				intVec2D{ intVec{ c }, intVec{ a,b,d,e } },
				intVec2D{ intVec{ d }, intVec{ a,b,c,e } },
				intVec2D{ intVec{ e }, intVec{ a,b,c,d } },
				intVec2D{ intVec{ a,b }, intVec{ c,d,e } },// 2+3
				intVec2D{ intVec{ a,c }, intVec{ b,d,e } },
				intVec2D{ intVec{ a,d }, intVec{ b,c,e } },
				intVec2D{ intVec{ a,e }, intVec{ b,c,d } },
				intVec2D{ intVec{ b,c }, intVec{ a,d,e } },
				intVec2D{ intVec{ b,d }, intVec{ a,c,e } },
				intVec2D{ intVec{ b,e }, intVec{ a,c,d } },
				intVec2D{ intVec{ c,d }, intVec{ a,b,e } },
				intVec2D{ intVec{ c,e }, intVec{ a,b,d } },
				intVec2D{ intVec{ d,e }, intVec{ a,b,c } },
				intVec2D{ intVec{ a }, intVec{ b }, intVec{ c,d,e } },// 1+1+3
				intVec2D{ intVec{ a }, intVec{ c }, intVec{ b,d,e } },
				intVec2D{ intVec{ a }, intVec{ d }, intVec{ b,c,e } },
				intVec2D{ intVec{ a }, intVec{ e }, intVec{ b,c,d } },
				intVec2D{ intVec{ b }, intVec{ c }, intVec{ a,d,e } },
				intVec2D{ intVec{ b }, intVec{ d }, intVec{ a,c,e } },
				intVec2D{ intVec{ b }, intVec{ e }, intVec{ a,c,d } },
				intVec2D{ intVec{ c }, intVec{ d }, intVec{ a,b,e } },
				intVec2D{ intVec{ c }, intVec{ e }, intVec{ a,b,d } },
				intVec2D{ intVec{ d }, intVec{ e }, intVec{ a,b,c } },
				intVec2D{ intVec{ e }, intVec{ a,b }, intVec{ c,d } },// 1+2+2
				intVec2D{ intVec{ d }, intVec{ a,b }, intVec{ c,e } },
				intVec2D{ intVec{ c }, intVec{ a,b }, intVec{ d,e } },
				intVec2D{ intVec{ e }, intVec{ a,c }, intVec{ b,d } },
				intVec2D{ intVec{ d }, intVec{ a,c }, intVec{ b,e } },
				intVec2D{ intVec{ b }, intVec{ a,c }, intVec{ d,e } },
				intVec2D{ intVec{ e }, intVec{ a,d }, intVec{ b,c } },
				intVec2D{ intVec{ c }, intVec{ a,d }, intVec{ b,e } },
				intVec2D{ intVec{ b }, intVec{ a,d }, intVec{ c,e } },
				intVec2D{ intVec{ d }, intVec{ a,e }, intVec{ b,c } },
				intVec2D{ intVec{ c }, intVec{ a,e }, intVec{ b,d } },
				intVec2D{ intVec{ b }, intVec{ a,e }, intVec{ c,d } },
				intVec2D{ intVec{ a }, intVec{ b,c }, intVec{ d,e } },
				intVec2D{ intVec{ a }, intVec{ b,d }, intVec{ c,e } },
				intVec2D{ intVec{ a }, intVec{ b,e }, intVec{ c,d } },
				intVec2D{ intVec{ c }, intVec{ d }, intVec{ e }, intVec{ a,b } },// 1+1+1+2
				intVec2D{ intVec{ b }, intVec{ d }, intVec{ e }, intVec{ a,c } },
				intVec2D{ intVec{ b }, intVec{ c }, intVec{ e }, intVec{ a,d } },
				intVec2D{ intVec{ b }, intVec{ c }, intVec{ d }, intVec{ a,e } },
				intVec2D{ intVec{ a }, intVec{ d }, intVec{ e }, intVec{ b,c } },
				intVec2D{ intVec{ a }, intVec{ c }, intVec{ e }, intVec{ b,d } },
				intVec2D{ intVec{ a }, intVec{ c }, intVec{ d }, intVec{ b,e } },
				intVec2D{ intVec{ a }, intVec{ b }, intVec{ e }, intVec{ c,d } },
				intVec2D{ intVec{ a }, intVec{ b }, intVec{ d }, intVec{ c,e } },
				intVec2D{ intVec{ a }, intVec{ b }, intVec{ c }, intVec{ d,e } },
				intVec2D{ intVec{ a }, intVec{ b }, intVec{ c }, intVec{ d }, intVec{ e } }// 1+1+1+1+1
		};
	case 1:
		return std::vector<intVec2D>{};
	default:
		std::cout << "Error: >5 bonds linked to an atom." << std::endl;
		throw "Error: >5 bonds linked to an atom."; // ------ throw error
	}
}
/* =========================================== */


/* =========================================== */
// important function: filter a molecule (or residue) to dup + res
// and give a new histogram to the next level
// note: the last 5 variable is transferred back
void filter(int Ntry, std::vector<strType> x0, std::vector<size_t> y0,
	std::unordered_map<strType, intVec2D> xID0, // input parameters
	std::vector<strType>& x, std::vector<size_t>& y,
	std::unordered_map<strType, intVec2D>& xID,
	intVec& dupN, std::vector<intVec2D>& dupIDs) {

	// generate the random pool used to pick a random fragment
	std::vector<size_t> randPool{};
	for (size_t i = 0; i < y0.size(); ++i) {
		for (size_t j = 0; j < y0[i]; ++j) {
			randPool.push_back(i);
		}
	}
	size_t PoolSize = randPool.size();
	// xNbond record #bonds for each histogram element, for speed
	std::vector<size_t> xNbond{};
	for (size_t i = 0; i < x0.size(); ++i) {
		xNbond.push_back(xID0[x0[i]][0].size());
	}

	intVec2D bondsTaken{}; // record which bonds have been filled
	// random pick fragments to add, until all empty slots gone
	int nUnchanged = 0;
	while (nUnchanged < Ntry) {
		nUnchanged++;
		// ------ 1 fragment process ------
		// random pick 1 fragment from the histogram, to add on

		intVec2D& IDs = xID0[x0[randPool[rand() % PoolSize]]];
		// try to add this chosen fragment to empty slot
		size_t istart = rand() % IDs.size();
		auto it = IDs.begin(); // each bonds-id associated with this inchi
		std::advance(it, istart);
		auto itStart = it;
		bool newFragAttached = false;
		while (true) {
			intVec linkedTo((*it).size(), -1);
			int whichCase = tryGunshot(*it, bondsTaken, linkedTo);
			if (whichCase == 2) {// totally covered
			}
			else if (whichCase == 1) { // totally available
				bondsTaken.push_back(*it);
				nUnchanged = 0;
				break;
			}
			else { // whichCase == 3
				intVec newFrag{}, alreadyTaken{};
				for (size_t i = 0; i < linkedTo.size(); ++i) {
					if (linkedTo[i] == -1) {
						newFrag.push_back((*it)[i]);
					}
					else {
						if (std::find(alreadyTaken.begin(), alreadyTaken.end(), linkedTo[i]) == alreadyTaken.end()) {
							alreadyTaken.push_back(linkedTo[i]);
							newFrag.insert(newFrag.end(), bondsTaken[linkedTo[i]].begin(), bondsTaken[linkedTo[i]].end());
						}
					}
				}
				size_t newFragNbond = newFrag.size();
				std::sort(newFrag.begin(), newFrag.end());
				for (size_t i = 0; i < xNbond.size(); ++i) {
					if (xNbond[i] == newFragNbond) {
						if (!isVecInVec2D(newFrag, xID0[x0[i]])) {
							continue; // didn't find
						}
						else { // found
							// delete already-existing overlapped from bondsTaken
							std::sort(linkedTo.begin(), linkedTo.end(), std::greater<>());
							for (size_t j = 0; j < linkedTo.size() - 1; ++j) {
								if (linkedTo[j] != -1
									&& linkedTo[j] != linkedTo[j + 1]) {
									bondsTaken.erase(bondsTaken.begin() + linkedTo[j]);
								}
							}
							if (linkedTo.back() != -1) {
								bondsTaken.erase(bondsTaken.begin() + linkedTo.back());
							}
							// add the new fragment
							bondsTaken.push_back(newFrag);
							newFragAttached = true;
							nUnchanged = 0;
							break;
						}
					}
				}
			}
			if (newFragAttached)
				break; // exist while
			it++;
			if (it == IDs.end())
				it = IDs.begin();
			if (it == itStart)
				break;
		}
		// ------ 1 fragment process ends ------
	}

	dupN.clear();
	dupIDs.clear();
	intVec2D res{};
	intVec resXid{}, dupXid{};
	for (intVec& frag : bondsTaken) {
		for (int i = 0; i < x0.size(); ++i) {
			if (isVecInVec2D(frag, xID0[x0[i]])) {
				// means found this fragment in "xIDthis"
				// index i uniquely represent a fragment

				if (std::find(resXid.begin(), resXid.end(), i) == resXid.end()) {
					resXid.push_back(i);
					res.push_back(frag);
				}
				else {
					auto dupXidFound = std::find(dupXid.begin(), dupXid.end(), i);
					if (dupXidFound == dupXid.end()) {
						dupXid.push_back(i);
						dupN.push_back(1);
						dupIDs.push_back(intVec2D{ frag });
					}
					else {
						auto idx = std::distance(dupXid.begin(), dupXidFound);
						dupN[idx]++;
						dupIDs[idx].push_back(frag);
					}
				}
				break;
			}
		}
	}


	//make the new histogram for the residue
	//new histogram consists of Rx, Ry, RxID
	std::vector<strType> RxPre{};
	std::vector<size_t> RyPre{};
	std::unordered_map<strType, intVec2D> RxIDPre{};
	for (size_t i = 0; i < res.size(); ++i) {

		//check every entry in x0 (namely fraghist.x)
		for (size_t j = 0; j < xNbond.size(); ++j) {
			// check if any item in fraghist.x is embraced by thisRes
			if (xNbond[j] <= res[i].size()) {
				intVec2D newIDset{};
				for (intVec ID : xID0[x0[j]]) {
					if (isEmbracing(res[i], ID)) {
						newIDset.push_back(ID);
					}
				}
				if (!newIDset.empty()) { // if smaller fragment is found
					auto RxFound = std::find(RxPre.begin(), RxPre.end(), x0[j]);
					if (RxFound == RxPre.end()) {
						// new inchi
						RxPre.push_back(x0[j]);
						RyPre.push_back(newIDset.size());
						RxIDPre[x0[j]] = newIDset;
					}
					else { // inchi already exists
						RyPre[std::distance(RxPre.begin(), RxFound)] += newIDset.size();
						intVec2D& it = RxIDPre[x0[j]];
						it.insert(it.end(), newIDset.begin(), newIDset.end());
					}
				}
			}
		}
	}
	// filter, only left with Ry >= 2
	x.clear();
	y.clear();
	xID.clear();
	for (size_t i = 0; i < RyPre.size(); ++i) {
		if (RyPre[i] >= 2) {
			x.push_back(RxPre[i]);
			y.push_back(RyPre[i]);
			xID[RxPre[i]] = RxIDPre[RxPre[i]];
		}
	}
}



/* =========================================== */
// construct a MOL instance from 
// (intVec& AtomsList, intVec2D& BondsList, int nfrag0)
MOL_BOND::MOL_BOND(intVec& AtomsList, intVec2D& BondsList, int nfrag0) {
	for (int i = 0; i < AtomsList.size(); ++i)
		Atoms[i] = ATOM{ AtomsList[i] };

	int a1, a2;
	for (int i = 0; i < BondsList.size(); ++i) {
		a1 = BondsList[i][0];
		a2 = BondsList[i][1];
		Bonds[i] = BOND{ BondsList[i][2],a1,a2 };
		Atoms[a1].BondsLinked.push_back(i);
		Atoms[a2].BondsLinked.push_back(i);
	}

	for (auto atom = Atoms.begin(); atom != Atoms.end(); atom++) {
		atom->second.toCut = getToCut(atom->second.BondsLinked);
		if (atom->second.BondsLinked.size() > 1) {
			AtomsHub.push_back(atom->first);
		}
	}
}


/* =========================================== */
// import data from external input file, transfering to class MOL_BOND
int MOL_BOND::importData(strType& molPATH, std::vector<MOL_BOND>& outputMolVec) {
	// read MolFile, and get rid of all H atoms
	MOB_IO_IreadMolFile::MolData molFileData 
		= MOB_IO_IreadMolFile::parse_mol_file(molPATH);
	if (!molFileData.good) {
		return -1; // file read fail
	}

	MOL_BOND mol0 = getMOL_from(molFileData, 1);
	// clean up: delete single atoms, single bonds,
	// and return a list of disconnected MOL_BOND.
	for (auto atom = mol0.Atoms.begin(); atom != mol0.Atoms.end(); ++atom) {
		if (atom->second.BondsLinked.size() == 0) {
			mol0.Atoms.erase(atom);
		}
	}
	std::vector<MOL_BOND> molList = mol0.split();
	for (MOL_BOND& mol : molList) {
		for (auto atom = mol.Atoms.begin(); atom != mol.Atoms.end(); atom++) {
			atom->second.toCut = getToCut(atom->second.BondsLinked);
			if (atom->second.BondsLinked.size() > 1) {
				mol.AtomsHub.push_back(atom->first);
			}
		}
	}


	if (molList.size() == 0) {
		return 0;
	}
	else if (molList.size() == 1) {
		molList[0].inchi = getInchi(molList[0]);
		molList[0].nAtomsIncludeH = getTotalNumAtoms(molList[0], molList[0].inchi);
		outputMolVec = molList;
		return 1;
	}
	else { // molList.size() >= 2
		for (MOL_BOND& mol : molList) {
			mol.inchi = getInchi(mol);
			mol.nAtomsIncludeH = getTotalNumAtoms(mol, mol.inchi);
		}
		outputMolVec = molList;
		return 2;
	}
}


/* =========================================== */
// cut and split, and return a vector of obtained MOL_BOND
std::vector<std::tuple<strType, strType, unsigned long long int>> MOL_BOND::cut_and_split() {
	MOL_BOND molcut = cutOnceWhole();
	std::vector<std::tuple<strType, strType, unsigned long long int>> histItem{};

	// split molcut, most of the following codes are directly copyed from
	// std::vector<MOL_BOND> MOL_BOND::split(), for sake of speed
	std::stack<int> sta{};
	int nextAtom, currentAtom, tBondID;
	intVec alreadyChecked{};
	bool found;

	MOL_BOND frag0{};
	unsigned long long int M0; // refer to math doc: weight value
	std::unordered_set<int> AtomsFrag0{}; // atoms (original order) included in fragment
	strType identifier{};
	intVec identifierVec{};

	auto it = molcut.Atoms.begin();
	while (it != molcut.Atoms.end()) {
		// --------------------- Atoms ---------------------
		found = find_x_in_V(it->first, alreadyChecked);
		if (found) { // if "it" is already checked, skip
			it++;
			continue;
		}
		else {
			// if "it" is a single bond, do not consider
			if (it->second.BondsLinked.size() <= 1) {
				tBondID = it->second.BondsLinked[0];
				nextAtom = molcut.Bonds[tBondID].Atom1id;
				if (nextAtom == it->first)
					nextAtom = molcut.Bonds[tBondID].Atom2id;
				if (molcut.Atoms[nextAtom].BondsLinked.size() <= 1) {
					if (molcut.Atoms[nextAtom].BondsLinked[0] == tBondID) {
						alreadyChecked.push_back(it->first);
						it++;
						continue;
					}
				}
			}
		}

		// initialise depth-first search
		frag0.Atoms.clear();
		frag0.Bonds.clear();
		AtomsFrag0.clear();
		identifierVec.clear();
		M0 = 1;
		frag0.Atoms[it->first] = it->second;
		for (int& i : it->second.BondsLinked) {
			nextAtom = molcut.Bonds[i].Atom1id;
			if (nextAtom == it->first)
				nextAtom = molcut.Bonds[i].Atom2id;
			sta.push(nextAtom);
		}
		alreadyChecked.push_back(it->first);

		// depth-first search
		while (!sta.empty()) {
			currentAtom = sta.top();
			sta.pop();
			alreadyChecked.push_back(currentAtom);
			frag0.Atoms[currentAtom] = molcut.Atoms[currentAtom];
			for (int& i : molcut.Atoms[currentAtom].BondsLinked) {
				nextAtom = molcut.Bonds[i].Atom1id;
				if (nextAtom == currentAtom)
					nextAtom = molcut.Bonds[i].Atom2id;
				if (!find_x_in_V(nextAtom, alreadyChecked))
					sta.push(nextAtom);
			}
		}

		// --------------------- Bonds ---------------------
		for (auto it2 = frag0.Atoms.begin(); it2 != frag0.Atoms.end(); ++it2) {
			for (int& j : it2->second.BondsLinked) {
				if (frag0.Bonds.find(j) == frag0.Bonds.end()) {
					// if j is not found, then add. This avoids repeats.
					BOND& jBond = molcut.Bonds[j];
					frag0.Bonds[j] = jBond;
					if (jBond.origAtom1id != -1) {
						AtomsFrag0.insert(jBond.origAtom1id);
						M0 *= getToCutInt(Atoms[jBond.origAtom1id].BondsLinked.size()
							- frag0.Atoms[jBond.Atom1id].BondsLinked.size());
					}
					if (jBond.origAtom2id != -1) {
						AtomsFrag0.insert(jBond.origAtom2id);
						M0 *= getToCutInt(Atoms[jBond.origAtom2id].BondsLinked.size()
							- frag0.Atoms[jBond.Atom2id].BondsLinked.size());
					}
					identifierVec.push_back(j);
				}
			}
		}
		identifier = getIdentifierStr(identifierVec);

		// calculate M0 further, contributed by atoms not inclued in this MOL
		for (int& i : AtomsHub) {
			if (AtomsFrag0.find(i) == AtomsFrag0.end()
				&& frag0.Atoms.find(i) == frag0.Atoms.end()) {
				M0 *= getToCutInt(Atoms[i].BondsLinked.size());
			}
		}

		histItem.push_back(std::make_tuple(getInchi(frag0), identifier, M0));
		it++;
	}
	return histItem;
}


/* =========================================== */
// calculate the total number of cutting schemes the original molecule has.
unsigned long long int MOL_BOND::calTotalM() {
	unsigned long long int M = 1;
	for (int& i : AtomsHub) {
		M *= getToCutInt(Atoms[i].BondsLinked.size());
	}
	return M;
}


/* =========================================== */
// main function in step 2: get 1 pathway, no matter short or not
PATHSINGLE MOL_BOND::get1pathway(FRAGHIST& fraghist, int Ntry_freeze) {
	PATHSINGLE path{};
	if (fraghist.y.size() == 0) {
		path.PAindex = static_cast<int>(Bonds.size()) - 1;
		return path;
	}

	int Ntry; // if Ntry times random events make no change, then stop
	if (Ntry_freeze < 0) {
		Ntry = static_cast<int>(Bonds.size()) * Ntry_freeze_amplifier;
	}
	else {
		Ntry = Ntry_freeze;
	}

	// get useful data from fraghist
	std::vector<size_t> y = fraghist.y;
	std::vector<strType> x = fraghist.x;
	std::unordered_map<strType, intVec2D> xID{};
	for (size_t i = 0; i < x.size(); ++i) {
		intVec2D xID_item{};
		for (strType idStr : fraghist.xID[x[i]]) {
			xID_item.push_back(retrieveBondsID(idStr));
		}
		xID[x[i]] = xID_item;
	}
	std::vector<strType> x0 = x;
	std::unordered_map<strType, intVec2D> xID0 = xID;
	intVec dupN{}; // one-to-one corresponding to index of fraghist.x
	std::vector<intVec2D> dupIDs{};

	// ===========================
	// ------ main function ------
	std::vector<strType> DupWrap{};
	intVec DupNWrap{};
	std::vector<intVec2D> DupIDsWrap{};
	while (true) {
		filter(Ntry, x, y, xID, // input parameters
			x, y, xID,
			dupN, dupIDs);
		for (size_t i = 0; i < dupN.size(); ++i) {
			size_t idx_x0 = findIDsIndexInX0(dupIDs[i][0], x0, xID0);
			auto it = std::find(DupWrap.begin(), DupWrap.end(), x0[idx_x0]);
			if (it == DupWrap.end()) {
				DupWrap.push_back(x0[idx_x0]);
				DupNWrap.push_back(dupN[i]);
				DupIDsWrap.push_back(dupIDs[i]);
			}
			else {
				auto idx = std::distance(DupWrap.begin(), it);
				DupNWrap[idx] += dupN[i];
				DupIDsWrap[idx].insert(DupIDsWrap[idx].end(),
					dupIDs[i].begin(), dupIDs[i].end());
			}
		}
		if (dupN.empty() || x.empty())
			break;
	}
	// ===========================
	// collect data
	// 1st: make the order vector: sort by DupNWrap (descending)
	std::vector<std::pair<size_t, int> > order(DupNWrap.size());
	size_t n = 0;
	for (auto it = DupNWrap.begin(); it != DupNWrap.end(); ++it, ++n)
		order[n] = std::make_pair(n, *it);
	std::sort(order.begin(), order.end(),
		[](auto& a, auto& b) { return a.second > b.second; });

	for (std::pair<size_t, int>& i : order) {
		path.dup.push_back(DupWrap[i.first]);
		path.dupN.push_back(DupNWrap[i.first]);
		path.dupIDs.push_back(DupIDsWrap[i.first]);
	}

	size_t paindex = Bonds.size() - 1;
	for (size_t i = 0; i < DupIDsWrap.size(); ++i) {
		paindex -= DupNWrap[i] * (DupIDsWrap[i][0].size() - 1);
	}
	path.PAindex = static_cast<int>(paindex);
	return path;
}


/* =========================================== */
// fint the minimal pathways (the final step)
void MOL_BOND::getMinPathways(PATHWAY_ASSEMBLY& PA,
	FRAGHIST& fraghist, int NtryPath, int Ntry_freeze) {
	PA.trivial_index = static_cast<int>(Bonds.size()) - 1;

	for (int i = 0; i < NtryPath; ++i) {
		PATHSINGLE path = get1pathway(fraghist, Ntry_freeze);
		if (path.PAindex < PA.index) {
			PA.paths.clear();
			if (!path.dup.empty()) {
				PA.paths.push_back(path);
			}
			PA.index = path.PAindex;
		}
		else if (path.PAindex == PA.index) {
			if (!path.dup.empty()) {
				bool isPathNew = true;
				for (size_t j = 0; j < PA.paths.size(); ++j) {
					if (isPathEqual(path, PA.paths[j])) {
						isPathNew = false;
						break;
					}
				}
				if (isPathNew) { // if path is new, then add it to PA.pahts
					PA.paths.push_back(path);
				}
			}
		}
	}
}



/* =========================================== */
//             private functions
/* =========================================== */
// randomly cut the molecule, return a MOL (original mol unchanged)
// only update .Atoms and .Bonds, thus not a proper MOL, only for "cut_and_split"
MOL_BOND MOL_BOND::cutOnceWhole() {
	MOL_BOND mol = MOL_BOND(*this);
	int newAtomID;
	ATOM newAtom;
	intVec2D currentCut;

	for (int& i : AtomsHub) {// go through every atoms that have more than 1 bond
		currentCut = Atoms[i].toCut[rand() % Atoms[i].toCut.size()]; // random choose cutting scheme for an atom

		if (currentCut.size() > 1) { // when there are actual cuts
			// go through this cutting scheme
			for (int j = 0; j < currentCut.size() - 1; j++) {// last cut should remain unchanged
				// add a new atom
				newAtomID = static_cast<int>(mol.Atoms.size());
				newAtom = ATOM();
				newAtom.AtomicNum = Atoms[i].AtomicNum;
				newAtom.BondsLinked = currentCut[j];
				mol.Atoms[newAtomID] = newAtom;

				// change bonds' end to this new atom
				for (int& k : currentCut[j]) {
					if (mol.Bonds[k].Atom1id == i) {
						mol.Bonds[k].Atom1id = newAtomID;
						mol.Bonds[k].origAtom1id = i;
					}
					else {
						mol.Bonds[k].Atom2id = newAtomID;
						mol.Bonds[k].origAtom2id = i;
					}
				}
			}
			for (int& k : currentCut.back()) {// indicate cut occurred, but atom didn't change
				if (mol.Bonds[k].Atom1id == i) {
					mol.Bonds[k].origAtom1id = i;
				}
				else {
					mol.Bonds[k].origAtom2id = i;
				}
			}

			// delete bonds-being-cut for this atom
			mol.Atoms[i].BondsLinked = currentCut.back();
		}
	}
	return mol;
}


/* =========================================== */
// splits the disconnected parts of mol, returns a vector of MOL 
// where each one of them is a disconnected part.
// Note that single bonds are neglected.
// only .Atoms and .Bonds updated
std::vector<MOL_BOND> MOL_BOND::split() {
	std::stack<int> sta{};
	int nextAtom, currentAtom, tBondID;
	std::vector<MOL_BOND> fragList{};
	MOL_BOND frag0{};
	intVec alreadyChecked{};
	bool found;

	auto it = Atoms.begin();
	while (it != Atoms.end()) {
		// --------------------- Atoms ---------------------
		found = find_x_in_V(it->first, alreadyChecked);
		if (found) { // if "it" is already checked, skip
			it++;
			continue;
		}
		else {
			// if "it" is a single bond, do not consider
			if (it->second.BondsLinked.size() <= 1) {
				tBondID = it->second.BondsLinked[0];
				nextAtom = Bonds[tBondID].Atom1id;
				if (nextAtom == it->first)
					nextAtom = Bonds[tBondID].Atom2id;
				if (Atoms[nextAtom].BondsLinked.size() <= 1) {
					if (Atoms[nextAtom].BondsLinked[0] == tBondID) {
						alreadyChecked.push_back(it->first);
						it++;
						continue;
					}
				}
			}
		}

		// initialise depth-first search
		frag0.Atoms.clear();
		frag0.Bonds.clear();
		frag0.Atoms[it->first] = it->second;
		for (int& i : it->second.BondsLinked) {
			nextAtom = Bonds[i].Atom1id;
			if (nextAtom == it->first)
				nextAtom = Bonds[i].Atom2id;
			sta.push(nextAtom);
		}
		alreadyChecked.push_back(it->first);

		// depth-first search
		while (!sta.empty()) {
			currentAtom = sta.top();
			sta.pop();
			alreadyChecked.push_back(currentAtom);
			frag0.Atoms[currentAtom] = Atoms[currentAtom];
			for (int& i : Atoms[currentAtom].BondsLinked) {
				nextAtom = Bonds[i].Atom1id;
				if (nextAtom == currentAtom)
					nextAtom = Bonds[i].Atom2id;
				if (!find_x_in_V(nextAtom, alreadyChecked))
					sta.push(nextAtom);
			}
		}

		// --------------------- Bonds ---------------------
		for (auto it2 = frag0.Atoms.begin(); it2 != frag0.Atoms.end(); ++it2) {
			for (int& j : it2->second.BondsLinked) {
				if (frag0.Bonds.find(j) == frag0.Bonds.end())
					// if j is not found, then add. This avoids repeats.
					frag0.Bonds[j] = Bonds[j];
			}
		}

		fragList.push_back(frag0);
		it++;
	}
	return fragList;
}
