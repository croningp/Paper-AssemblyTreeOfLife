#include "MOB_IO_O.h"
#include "MOB_IO_IconvertToMOB.h"
#include <iostream>
#include <fstream>


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


/* =========================================== */
// The final output to console.
void outputConsole(MOL_BOND& mol, PATHWAY_ASSEMBLY& PA,
	strType& outPATH, strType& outPATH_hist,
	long long int Nstep1, long long int Nstep2, 
	long long int nTotalTry, double ElapsedTime) {
	std::cout << "==================================================" << std::endl;
	std::cout << "============ Printing step 2 results =============" << std::endl;
	std::cout << "==================================================" << std::endl;
	std::cout << "The molecule to be analysed is:" << std::endl
		<< mol.inchi << std::endl;
	std::cout << "Output stored in: " << outPATH << " and " << outPATH_hist << std::endl;
	std::cout << ". Parameter Nstep1 = " << Nstep1 << std::endl;
	std::cout << ". Parameter Nstep2 = " << Nstep2 << std::endl;
	std::cout << ". " << nTotalTry << " pathways tried." << std::endl;
	std::cout << ". Elapsed time: " << ElapsedTime << " s" << std::endl;
	std::cout << std::endl;
	std::cout << "Number of bonds (excluding H): " << mol.Bonds.size() << std::endl;
	std::cout << "Number of pathways: " << PA.paths.size() << std::endl;
	std::cout << "Assembly index is: " << PA.index << std::endl;


	int i = 0;
	for (PATHSINGLE& path : PA.paths) {
		i++;
		std::cout << "=================" << std::endl;
		std::cout << "=== Pathway " << i << " ===" << std::endl;
		std::cout << "=================" << std::endl;

		for (strType& inchi : path.dup) {
			std::cout << inchi << std::endl;
		}
		std::cout << std::endl;

		std::cout << "How many times each duplicates:" << std::endl;
		for (size_t i = 0; i < path.dup.size(); ++i) {
			std::cout << path.dupN[i] << " :  ";
			auto& vec2 = path.dupIDs[i];
			for (auto& vec : vec2) {
				std::cout << "(";
				for (size_t j = 0; j < vec.size() - 1; ++j) {
					std::cout << vec[j] << ",";
				}
				std::cout << vec.back() << ")  ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl << std::endl;
	}
}


/* =========================================== */
// output the fragment histogram
void outputFile_hist(MOL_BOND& mol, FRAGHIST& fragHist, 
	strType& outPATH_hist,
	long long int Nstep1, double ElapsedTime) {

	std::ofstream pathVecFile;
	pathVecFile.open(outPATH_hist);
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
	pathVecFile << "==================================================" << std::endl;
	pathVecFile << "========= Printing fragments histogram ===========" << std::endl;
	pathVecFile << "==================================================" << std::endl;
	pathVecFile << "Parameter Nstep1 = " << Nstep1 << std::endl;
	pathVecFile << "Elapsed time: " << ElapsedTime << " s" << std::endl << std::endl;

	int nGroup = 5;
	for (size_t i = 0; i < fragHist.x.size(); i += nGroup) {
		for (size_t j = 0; j < nGroup; ++j) {
			if (i+j < fragHist.x.size()) {
				pathVecFile << fragHist.x[i + j] << std::endl;
			}
		}
		for (size_t j = 0; j < nGroup; ++j) {
			if (i + j < fragHist.x.size()) {
				pathVecFile << fragHist.y[i + j] << std::endl;
			}
		}
		pathVecFile << std::endl;

	}
	pathVecFile.close();
}


/* =========================================== */
// Print instruction in console
void printInstruction() {
	strType exeName = "[this-exe-file-name].exe";
	std::cout << "How to run the exe-file in command line:" << std::endl
		<< "1. Enter \'" << exeName << "\'," << std::endl
		<< "     it will analyse the mol-file \'AnalyseMeDefault.mol\';" << std::endl
		<< "2. Enter \'" << exeName << " [mol-file-name].mol\'," << std::endl
		<< "     it will analyse the mol-file \'[mol-file-name].mol\'." << std::endl
		<< "Parameters: "
		<< "There are 2 parameters that you can input in order either 0, 1 or 2 of them." << std::endl
		<< " .1st parameter Nstep2 means how many pathways to try to find the minimal" << std::endl
		<< "     (-1 is defualt, meaning never stop)" << std::endl
		<< " .2nd parameter Nstep1 means how many cutting schemes to try to obtain fragment histogram" << std::endl
		<< "     (-1 is defualt, meaning max(" << Nstep1_Min << ", 1% of possible schemes))" << std::endl << std::endl;
}
