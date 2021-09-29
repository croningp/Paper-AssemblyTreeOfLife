/* The main function.
   It calculates the assembly space for one single molecule.
*/

#include "MOB.h"
#include "FRAGHIST.h"
#include "MOB_IO_O.h"
#include <iostream>
#include <ctime>
#include <set>
#include <chrono>
typedef MOL_BOND OBJ;
void printInstruction(); // Print instruction in console


int main(int argc, char* argv[]) {
	srand(static_cast<unsigned int>(time(0)));
	//srand(99); // fix random seed ----- for test -----

	std::cout << std::endl << versionStr << std::endl << std::endl;
	printInstruction();


	// ! Command line wrong input: too many arguments
	if (argc > 4) {
		std::cout << std::endl << "! Command line wrong input." << std::endl;
		return 0;
	}


	long long int Nstep2 = -1, Nstep1 = -1;
	strType molPATH{}, outPATH{}, outPATH_hist{};
	switch (argc)
	{
	case 1:
		molPATH = "AnalyseMeDefault.mol";
		outPATH_hist = "AnalyseMeDefault_histogram.txt";
		outPATH = "AnalyseMeDefault_pathway.txt";
		break;
	case 2:
		molPATH = argv[1];
		outPATH_hist = molPATH.substr(0, molPATH.size() - 4) + "_histogram.txt";
		outPATH = molPATH.substr(0, molPATH.size() - 4) + "_pathway.txt";
		break;
	case 3:
		molPATH = argv[1];
		outPATH_hist = molPATH.substr(0, molPATH.size() - 4) + "_histogram.txt";
		outPATH = molPATH.substr(0, molPATH.size() - 4) + "_pathway.txt";
		try {
			Nstep2 = std::stoll(argv[2]);
		}
		catch (...) {
			std::cout << "! Parameters wrong input." << std::endl;
			return 0;
		}
		break;
	case 4:
		molPATH = argv[1];
		outPATH_hist = molPATH.substr(0, molPATH.size() - 4) + "_histogram.txt";
		outPATH = molPATH.substr(0, molPATH.size() - 4) + "_pathway.txt";
		try {
			Nstep2 = std::stoll(argv[2]);
			Nstep1 = std::stoll(argv[3]);
		}
		catch (...) {
			std::cout << "! Parameters wrong input." << std::endl;
			return 0;
		}
	}


	PATHWAY_ASSEMBLY PA{};
	// import data from mol file
	OBJ obj{};
	std::vector<OBJ> objList{};
	int intMsg = obj.importData(molPATH, objList);
	switch (intMsg) {
	case -1:
		std::cout << "! Cannot read mol-file." << std::endl;
		if (argc == 1) {
			system("pause");
		}
		return 0;
	case 0:
		std::cout << Err_Msg_900 << std::endl;
		outputFile(obj, PA, 900, outPATH, outPATH_hist, 0, 0, 0, 0);
		if (argc == 1) {
			system("pause");
		}
		return 0;
	case 2:
		std::cout << Err_Msg_902 << std::endl;
		outputFile(obj, PA, 902, outPATH, outPATH_hist, 0, 0, 0, 0);
		if (argc == 1) {
			system("pause");
		}
		return 0;
	case 1:
		// mol-file read & #mol = 1. Good input, can proceed.
		break;
	}
	obj = objList[0];

	unsigned long long int Mt = obj.calTotalM(); // total number of possible cutting schemes
	if (Nstep1 == -1) {
		Nstep1 = static_cast<long long int> (Mt / 100);
		if (Nstep1 < Nstep1_Min) {
			Nstep1 = Nstep1_Min;
		}
	}
	outputFile(obj, PA, 1, outPATH, outPATH_hist, Nstep1, Nstep2, 0, 0);


	std::chrono::time_point<std::chrono::steady_clock> startT, step1T, endT;
	startT = std::chrono::steady_clock::now();
	std::chrono::duration<double> running_time;
	double ElapsedTime;
	std::cout << std::endl << "Please wait..." << std::endl 
		<< " ." << Nstep1 << " steps to try in step 1..." << std::endl;
	step1T = std::chrono::steady_clock::now();
	// Monte Carlo step 1: generate fragments histogram
	std::vector<std::tuple<strType, strType, unsigned long long int>> histItem{};
	FRAGHIST fragHist{};
	for (long long int i = 0; i < Nstep1; ++i) {
		histItem = obj.cut_and_split();
		fragHist.addOnce(histItem);
		endT = std::chrono::steady_clock::now();
		running_time = endT - step1T;
		if (running_time.count() > 5) {
			std::cout << " >" << i;
			step1T = endT;
		}
	}
	fragHist.getHistOverall(Mt);
	fragHist.getFinalHist(true, Hist_Err_Threshold);
	endT = std::chrono::steady_clock::now();
	running_time = endT - startT;
	ElapsedTime = running_time.count();
	outputFile_hist(obj, fragHist, outPATH_hist, Nstep1, ElapsedTime);
	std::cout << std::endl
		<< "=============================" << std::endl
		<< "====== Step 1 is Done. ======" << std::endl
		<< "=============================" << std::endl;


	// ---------------------------------------------------------------
	// Monte Carlo step 2
	std::cout << ">";
	long long int Nloop = 0;
	if (Nstep2 == -1) {
		Nloop = -2;
	}

	unsigned long long int nTotalTry = 0;
	while (Nloop < Nstep2) {
		nTotalTry += Nstep2_SaveWindow;
		Nloop += Nstep2_SaveWindow;
		for (int i = 0; i < int(Nstep2_SaveWindow/100); ++i) {
			obj.getMinPathways(PA, fragHist, 100);
			std::cout << ">";
		}
		endT = std::chrono::steady_clock::now();
		running_time = endT - startT;
		ElapsedTime = running_time.count();
		// Write final output
		outputFile(obj, PA, 2, outPATH, outPATH_hist,
			Nstep1, Nstep2, nTotalTry, ElapsedTime);
		std::cout << std::endl << "New results exported, " << nTotalTry << " pathways tried." << std::endl;
		outputConsole(obj, PA, outPATH, outPATH_hist,
			Nstep1, Nstep2, nTotalTry, ElapsedTime);

		if (Nstep2 == -1) {
			Nloop = -2;
		}
	}
	// ---------------------------------------------------------------
	// ---------------------------------------------------------------
	return 1;
}
