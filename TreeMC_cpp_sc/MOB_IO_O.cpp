#include "MOB_IO_O.h"
#include "MOB_IO_IconvertToMOB.h"
#include <algorithm>


// 1st stage of output, just write down to files that the program is running.
void outputFile_stage1(MOL_BOND& mol, strType& outPATH,
	long long int Nstep1, long long int Nstep2) {
	std::ofstream pathVecFile;
	pathVecFile.open(outPATH);
	pathVecFile << versionStr << std::endl << std::endl;
	pathVecFile << "The molecule to be analysed is:" << std::endl
		<< mol.inchi << std::endl;

	// write the indexing of atoms and bonds of this MOL
	intVec temp;
	for (auto it = mol.Atoms.begin(); it != mol.Atoms.end(); ++it) {
		pathVecFile << "    ||  " << getAtomicSym2Char(it->second.AtomicNum) << "  "
			<< "#" << it->first << "  -->" << "  ";
		temp = it->second.BondsLinked;
		for (int i = 0; i < temp.size() - 1; ++i) {
			pathVecFile << temp[i] << ", ";
		}
		pathVecFile << temp.back() << std::endl;
	}
	for (auto it = mol.Bonds.begin(); it != mol.Bonds.end(); ++it) {
		pathVecFile << "    ||      " << it->first << "  "
			<< "*" << it->second.BondType << " :  "
			<< "#" << it->second.Atom1id << " -- "
			<< "#" << it->second.Atom2id << std::endl;
	}


	pathVecFile << std::endl << ". Total number of atoms: "
		<< mol.nAtomsIncludeH << std::endl;
	pathVecFile << ". When hydrogen atoms (H) are excluded:" << std::endl;
	pathVecFile << "Number of atoms: " << mol.Atoms.size() << std::endl;
	pathVecFile << "Number of bonds: " << mol.Bonds.size() << std::endl;
	pathVecFile << std::endl;
	pathVecFile << ". Parameter Nstep1 = " << Nstep1 << std::endl;
	pathVecFile << ". Parameter Nstep2 = " << Nstep2 << std::endl;

	pathVecFile << std::endl << "Please wait until it generates results..."
		<< std::endl;
	pathVecFile.close();
}



/* =========================================== */
// The final output to files.
int outputFile(MOL_BOND& mol, PATHWAY_ASSEMBLY& PA, int stage,
	strType& outPATH, strType& outPATH_hist,
	long long int Nstep1, long long int Nstep2,
	long long int nTotalTry, double ElapsedTime) {
	std::ofstream pathVecFile;

	// ! Incorrect input: mol-file contains no valid molecule
	if (stage == 900) {
		pathVecFile.open(outPATH);
		pathVecFile << Err_Msg_900 << std::endl;
		return 1;
	}

	// ! Incorrect input: the number of molecules in the mol-file >= 2
	if (stage == 902) {
		pathVecFile.open(outPATH);
		pathVecFile << Err_Msg_902 << std::endl;
		return 1;
	}


	// output stage 1, just write the info of this molecule to the output file.
	if (stage == 1) {
		outputFile_stage1(mol, outPATH, Nstep1, Nstep2);
		return 1;
	}


	// 2nd stage of output, write down final output to files.
	if (stage == 2) {
		outputFile_stage1(mol, outPATH, Nstep1, Nstep2);
		pathVecFile.open(outPATH, std::fstream::app);
		pathVecFile << "==================================================" << std::endl;
		pathVecFile << "============ Printing step 2 results =============" << std::endl;
		pathVecFile << "==================================================" << std::endl;
		pathVecFile << ". " << nTotalTry << " pathways tried." << std::endl;
		pathVecFile << ". Elapsed time: " << ElapsedTime << " s" << std::endl;
		pathVecFile << std::endl;
		pathVecFile << "Number of pathways: " << PA.paths.size() << std::endl;
		pathVecFile << "Assembly index is: " << PA.index << std::endl << std::endl;


		int i = 0;
		for (PATHSINGLE& path : PA.paths) {
			i++;
			pathVecFile << "=================" << std::endl;
			pathVecFile << "=== Pathway " << i << " ===" << std::endl;
			pathVecFile << "=================" << std::endl;

			for (strType& inchi : path.dup) {
				pathVecFile << inchi << std::endl;
			}
			pathVecFile << std::endl;

			pathVecFile << "How many times each duplicates:" << std::endl;
			for (size_t i = 0; i < path.dup.size(); ++i) {
				pathVecFile << path.dupN[i] << " :  ";
				auto& vec2 = path.dupIDs[i];
				for (auto& vec : vec2) {
					pathVecFile << "(";
					for (size_t j = 0; j < vec.size() - 1; ++j) {
						pathVecFile << vec[j] << ",";
					}
					pathVecFile << vec.back() << ")  ";
				}
				pathVecFile << std::endl;
			}
			pathVecFile << std::endl << std::endl;
		}
		pathVecFile.close();
		return 1;
	}
	return 1;
}


///* =========================================== */
//// The final output to console.
//void outputConsole(MOL_BOND& mol, PATHWAY_ASSEMBLY& PA,
//	strType& outPATH, strType& outPATH_hist,
//	long long int Nstep1, long long int Nstep2,
//	long long int nTotalTry, double ElapsedTime) {
//	std::cout << "==================================================" << std::endl;
//	std::cout << "============ Printing step 2 results =============" << std::endl;
//	std::cout << "==================================================" << std::endl;
//	std::cout << "The molecule to be analysed is:" << std::endl
//		<< mol.inchi << std::endl;
//	std::cout << "Output stored in: " << outPATH << " and " << outPATH_hist << std::endl;
//	std::cout << ". Parameter Nstep1 = " << Nstep1 << std::endl;
//	std::cout << ". Parameter Nstep2 = " << Nstep2 << std::endl;
//	std::cout << ". " << nTotalTry << " pathways tried." << std::endl;
//	std::cout << ". Elapsed time: " << ElapsedTime << " s" << std::endl;
//	std::cout << std::endl;
//	std::cout << "Number of bonds (excluding H): " << mol.Bonds.size() << std::endl;
//	std::cout << "Number of pathways: " << PA.paths.size() << std::endl;
//	std::cout << "Assembly index is: " << PA.index << std::endl;
//
//
//	int i = 0;
//	for (PATHSINGLE& path : PA.paths) {
//		i++;
//		std::cout << "=================" << std::endl;
//		std::cout << "=== Pathway " << i << " ===" << std::endl;
//		std::cout << "=================" << std::endl;
//
//		for (strType& inchi : path.dup) {
//			std::cout << inchi << std::endl;
//		}
//		std::cout << std::endl;
//
//		std::cout << "How many times each duplicates:" << std::endl;
//		for (size_t i = 0; i < path.dup.size(); ++i) {
//			std::cout << path.dupN[i] << " :  ";
//			auto& vec2 = path.dupIDs[i];
//			for (auto& vec : vec2) {
//				std::cout << "(";
//				for (size_t j = 0; j < vec.size() - 1; ++j) {
//					std::cout << vec[j] << ",";
//				}
//				std::cout << vec.back() << ")  ";
//			}
//			std::cout << std::endl;
//		}
//		std::cout << std::endl << std::endl;
//	}
//}
//
//
///* =========================================== */
//// output the fragment histogram
//void outputFile_hist(MOL_BOND& mol, FRAGHIST& fragHist,
//	strType& outPATH_hist,
//	long long int Nstep1, double ElapsedTime) {
//
//	std::ofstream pathVecFile;
//	pathVecFile.open(outPATH_hist);
//	pathVecFile << versionStr << std::endl << std::endl;
//	pathVecFile << "The molecule to be analysed is:" << std::endl
//		<< mol.inchi << std::endl;
//
//	// write the indexing of atoms and bonds of this MOL
//	intVec temp;
//	for (auto it = mol.Atoms.begin(); it != mol.Atoms.end(); ++it) {
//		pathVecFile << "    ||  " << getAtomicSym2Char(it->second.AtomicNum) << "  "
//			<< "#" << it->first << "  -->" << "  ";
//		temp = it->second.BondsLinked;
//		for (int i = 0; i < temp.size() - 1; ++i) {
//			pathVecFile << temp[i] << ", ";
//		}
//		pathVecFile << temp.back() << std::endl;
//	}
//	for (auto it = mol.Bonds.begin(); it != mol.Bonds.end(); ++it) {
//		pathVecFile << "    ||      " << it->first << "  "
//			<< "*" << it->second.BondType << " :  "
//			<< "#" << it->second.Atom1id << " -- "
//			<< "#" << it->second.Atom2id << std::endl;
//	}
//
//
//	pathVecFile << std::endl << ". Total number of atoms: "
//		<< mol.nAtomsIncludeH << std::endl;
//	pathVecFile << ". When hydrogen atoms (H) are excluded:" << std::endl;
//	pathVecFile << "Number of atoms: " << mol.Atoms.size() << std::endl;
//	pathVecFile << "Number of bonds: " << mol.Bonds.size() << std::endl;
//	pathVecFile << std::endl;
//	pathVecFile << "==================================================" << std::endl;
//	pathVecFile << "========= Printing fragments histogram ===========" << std::endl;
//	pathVecFile << "==================================================" << std::endl;
//	pathVecFile << "Parameter Nstep1 = " << Nstep1 << std::endl;
//	pathVecFile << "Elapsed time: " << ElapsedTime << " s" << std::endl << std::endl;
//
//	int nGroup = 5;
//	for (size_t i = 0; i < fragHist.x.size(); i += nGroup) {
//		for (size_t j = 0; j < nGroup; ++j) {
//			if (i + j < fragHist.x.size()) {
//				pathVecFile << fragHist.x[i + j] << std::endl;
//			}
//		}
//		for (size_t j = 0; j < nGroup; ++j) {
//			if (i + j < fragHist.x.size()) {
//				pathVecFile << fragHist.y[i + j] << std::endl;
//			}
//		}
//		pathVecFile << std::endl;
//
//	}
//	pathVecFile.close();
//}

void outputFile_histWhole(MOL_BOND& mol, FRAGHIST& fragHist,
	strType& outPATH_histWhole,
	long long int Nstep1, double ElapsedTime) {

	std::ofstream pathVecFile;
	pathVecFile.open(outPATH_histWhole);
	pathVecFile << versionStr << std::endl << std::endl;
	pathVecFile << "The molecule analysed is:" << std::endl
		<< mol.inchi << std::endl;

	// write the indexing of atoms and bonds of this MOL
	intVec temp;
	for (auto it = mol.Atoms.begin(); it != mol.Atoms.end(); ++it) {
		pathVecFile << "    ||  " << getAtomicSym2Char(it->second.AtomicNum) << "  "
			<< "#" << it->first << "  -->" << "  ";
		temp = it->second.BondsLinked;
		for (int i = 0; i < temp.size() - 1; ++i) {
			pathVecFile << temp[i] << ", ";
		}
		pathVecFile << temp.back() << std::endl;
	}
	for (auto it = mol.Bonds.begin(); it != mol.Bonds.end(); ++it) {
		pathVecFile << "    ||      " << it->first << "  "
			<< "*" << it->second.BondType << " :  "
			<< "#" << it->second.Atom1id << " -- "
			<< "#" << it->second.Atom2id << std::endl;
	}


	pathVecFile << std::endl << ". Total number of atoms: "
		<< mol.nAtomsIncludeH << std::endl;
	pathVecFile << ". When hydrogen atoms (H) are excluded:" << std::endl;
	pathVecFile << "Number of atoms: " << mol.Atoms.size() << std::endl;
	pathVecFile << "Number of bonds: " << mol.Bonds.size() << std::endl;
	pathVecFile << std::endl;
	pathVecFile << "Parameter Nstep1 = " << Nstep1 << std::endl;
	pathVecFile << "Elapsed time: " << ElapsedTime << " s" << std::endl;
	pathVecFile << "==================================================" << std::endl;
	pathVecFile << "======= The semi-whole fragments histogram =======" << std::endl;
	pathVecFile << "==================================================" << std::endl;
	pathVecFile << "$" << std::endl;

	for (size_t i = 0; i < fragHist.x_orig.size(); ++i) {
		pathVecFile << fragHist.x_orig[i] << std::endl;
		pathVecFile << std::round(fragHist.y_mean_orig[i]) << ":" << std::endl;
		std::set<strType> IDs = fragHist.xID_orig[fragHist.x_orig[i]];
		for (strType id0 : IDs) {
			pathVecFile << "(" << id0 << ")  ";
		}
		pathVecFile << std::endl << std::endl;
	}
	pathVecFile.close();
}


/* =========================================== */
// Print instruction in console
//void printInstruction() {
//	strType exeName = "[this-exe-file-name].exe";
//	std::cout << "How to run the exe-file in command line:" << std::endl
//		<< "1. Enter \'" << exeName << "\'," << std::endl
//		<< "     it will analyse the mol-file \'AnalyseMeDefault.mol\';" << std::endl
//		<< "2. Enter \'" << exeName << " [mol-file-name].mol\'," << std::endl
//		<< "     it will analyse the mol-file \'[mol-file-name].mol\'." << std::endl
//		<< "Parameters: "
//		<< "There are 2 parameters that you can input in order either 0, 1 or 2 of them." << std::endl
//		<< " .1st parameter Nstep2 means how many pathways to try to find the minimal" << std::endl
//		<< "     (-1 is defualt, meaning never stop)" << std::endl
//		<< " .2nd parameter Nstep1 means how many cutting schemes to try to obtain fragment histogram" << std::endl
//		<< "     (-1 is defualt, meaning max(" << Nstep1_Min << ", 1% of possible schemes))" << std::endl << std::endl;
//}


/* =========================================== */
void outputAllPaths(std::vector<PATHWAY_ASSEMBLY>& PAvec, strType FilePath) {
	std::ofstream f;
	f.open(FilePath);
	for (auto& PA : PAvec) {
		if (!PA.paths.empty()) {
			f << "====================" << std::endl;
			f << "Index = " << PA.paths[0].PAindex << std::endl;
			for (auto& path : PA.paths) {
				f << "$" << std::endl;
				for (size_t i = 0; i < path.dup.size(); ++i) {
					f << path.dup[i] << std::endl;
					f << path.dupN[i] << ", " 
						<< path.dupIDs[i][0].size() << std::endl;
				}
				f << std::endl;
			}
		}
	}
	f.close();
}

/* =========================================== */
void outputTree(std::vector<PATHWAY_ASSEMBLY>& PAvec,
	strType FilePath, int ThresholdPA, strType orderBy) {
	std::ofstream f;
	f.open(FilePath);

	if (ThresholdPA >= PAvec.size() or ThresholdPA < 0) {
		ThresholdPA = static_cast<int>(PAvec.size()) - 1;
	}
	std::unordered_map<strType, intVec> TreeRaw{};
	for (size_t i = 1; i <= ThresholdPA; ++i) {
		for (PATHSINGLE& path : PAvec[i].paths) {
			for (size_t j = 0; j < path.dup.size(); ++j) {
				auto it = TreeRaw.insert({ path.dup[j], intVec{1, 0} });
				if (!it.second) {
					(it.first->second[0])++;
				}
				else {
					it.first->second[1] = static_cast<int>(path.dupIDs[j][0].size());
				}
			}
		}
	}

	// order the inchis
	std::vector<strType> inchi{};
	intVec Nrepeat{}, Nbonds{};
	for (auto& it : TreeRaw) {
		inchi.push_back(it.first);
		Nrepeat.push_back(it.second[0]);
		Nbonds.push_back(it.second[1]);
	}
	std::vector<std::pair<size_t, int> > order(Nrepeat.size());
	size_t n = 0;
	if (orderBy == "size") {
		for (auto it = Nbonds.begin(); it != Nbonds.end(); ++it, ++n)
			order[n] = std::make_pair(n, *it);
	}
	else if (orderBy == "nRepeats") {
		for (auto it = Nrepeat.begin(); it != Nrepeat.end(); ++it, ++n)
			order[n] = std::make_pair(n, *it);
	}
	else {
		throw "OutputTree function wrong input.";
		std::cout << "OutputTree function wrong input." << std::endl;
		return;
	}
	std::sort(order.begin(), order.end(),
		[](auto& a, auto& b) { return a.second > b.second; });
	std::vector<strType> inchi0{};
	intVec Nrepeat0{}, Nbonds0{};
	for (std::pair<size_t, int>& i : order) {
		inchi0.push_back(inchi[i.first]);
		Nrepeat0.push_back(Nrepeat[i.first]);
		Nbonds0.push_back(Nbonds[i.first]);
	}

	for (size_t i = 0; i < inchi0.size(); ++i) {
		f << inchi0[i] << std::endl;
		f << Nrepeat0[i] << ", " << Nbonds0[i] << std::endl;
	}
	f.close();
}
