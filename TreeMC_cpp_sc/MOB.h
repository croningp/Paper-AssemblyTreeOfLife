/* It defines MOL class and associated functions.
   MOL includes class ATOM and class BOND.
*/
#pragma once
#include "MOB_Def.h"
#include "FRAGHIST.h"
#include <unordered_map>


/* =========================================== */
// to record one single pathway (not necessarily minimal)
class PATHSINGLE {
public:
	std::vector<strType> dup{}; // vector of inchis of duplicated fragments
	intVec dupN{}; // how many times each duplication duplicates
	std::vector<intVec2D> dupIDs{}; // the detailed identifier of each fragment
	int PAindex = -1; // assembly index for this pathway

	PATHSINGLE() {};
};
class PATHWAY_ASSEMBLY {
public:
	std::vector<PATHSINGLE> paths{};
	int index = INT_MAX; // assembly index for this molecule
	int trivial_index = -1; // mol.Bonds.size() - 1;

	PATHWAY_ASSEMBLY() {};
};


/* =========================================== */
struct ATOM {
	int AtomicNum = 0;
	intVec BondsLinked{};
	std::vector<std::vector<intVec>> toCut{}; // store cutSchemes for each atom
};
struct BOND {
	int BondType = 0, Atom1id = 0, Atom2id = 0;
	int origAtom1id = -1, origAtom2id = -1;
	// record the original id of the replaced atom
};


/* =========================================== */
class MOL_BOND {
public:
	std::unordered_map<int, ATOM> Atoms{};
	std::unordered_map<int, BOND> Bonds{};

	intVec AtomsHub{}; // atoms with >=2 bonds
	strType inchi{};
	int nAtomsIncludeH = -1;


	// ------ functions ------
	MOL_BOND() {};
	MOL_BOND(intVec&, intVec2D&, int);

	int importData(strType& molPATH, std::vector<MOL_BOND>& outputMolVec);
	std::vector<std::tuple<strType, strType, unsigned long long int>> cut_and_split();
	std::vector<std::tuple<strType, strType, unsigned long long int>> histItemOrigMol();
	unsigned long long int calTotalM();
	PATHSINGLE get1pathway(FRAGHIST& fraghist, int Nmol, int Ntry_freeze = -1);
	void getMinPathways(std::unordered_map<strType, intVec>& TreeRaw, 
		PATHWAY_ASSEMBLY& PA,
		FRAGHIST& fraghist, int Nmol, int NtryPath, int Ntry_freeze = -1);
	// -- PCTree associated function --
	int addMol(strType molFile);
	std::vector<PATHWAY_ASSEMBLY> collectPaths(FRAGHIST& fraghist, int Nmol, int NtryPath, int Ntry_freeze = -1);
	//

private:
	std::vector<MOL_BOND> split();
	MOL_BOND cutOnceWhole();
};


strType getIdentifierStr(intVec& identifierVec);
intVec retrieveBondsID(strType identifierStr);
bool isPathEqual(PATHSINGLE& path1, PATHSINGLE& path2);