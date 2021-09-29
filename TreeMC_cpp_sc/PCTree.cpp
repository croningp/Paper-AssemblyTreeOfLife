#include "MOB.h"
#include "FRAGHIST.h"
#include "FRAGHIST_EXTERNAL.h"
#include "MOB_IO_O.h"
#include <iostream>
#include <fstream>
#include <chrono>
typedef MOL_BOND OBJ;


int main(int argc, char* argv[]) {
	unsigned int tt = static_cast<unsigned int>(time(0));
	//std::cout << "random seed: " << tt << std::endl; ----- for test -----
	srand(tt);
	//srand(99); // fix random seed ----- for test ----- seed = 99

	strType ToDoFilePath;
	if (argc == 1) {
		ToDoFilePath = "ToDo.txt";
	}
	else if (argc == 2) {
		ToDoFilePath = argv[1];
	}
	else {
		std::cout << "! Wrong: too many arguments." << std::endl;
		return 0;
	}
	std::ifstream ToDoFile(ToDoFilePath);
	if (ToDoFile.fail()) {
		std::cout << "! Wrong: ToDo-file opening fails." << std::endl;
		if (argc == 1)
			system("pause");
		return 0;
	}
	strType ToDoFolderPath;
	size_t i = ToDoFilePath.find_last_of('\\');
	if (i == std::string::npos) {
		ToDoFolderPath = ".\\";
	}
	else {
		ToDoFolderPath = ToDoFilePath.substr(0, i+1);
	}

	int Nstep2 = -1;
	std::vector<long long int> Nstep1vec{};
	std::vector<strType> MolTodo{};
	strType line{};
	std::getline(ToDoFile, line);
	i = line.find("=");
	if (i != std::string::npos)
		Nstep2 = std::stoi(line.substr(i + 1));
	while (std::getline(ToDoFile, line)) {
		if (line != "$") {
			MolTodo.push_back(line);
			Nstep1vec.push_back(-1);
			std::getline(ToDoFile, line);
			i = line.find("=");
			if (i != std::string::npos)
				Nstep1vec.back() = std::stoll(line.substr(i + 1));
		}
		else {
			break;
		}
	}
	ToDoFile.close();

	if (!(Nstep2 == -1 || Nstep2 > 0)) {
		std::cout << "! Wrong: Nstep2 value" << std::endl;
		if (argc == 1)
			system("pause");
		return 0;
	}
	std::ifstream f;
	for (size_t j = 0; j < Nstep1vec.size(); ++j) {
		strType temp = ToDoFolderPath + MolTodo[j] + "_histWhole.txt";
		f.open(temp);
		if (Nstep1vec[j] == 0) { // use already-existing hist-file
			if (f.fail()) {
				std::cout << "! Wrong: " << temp << " does not exist." << std::endl;
				if (argc == 1)
					system("pause");
				return 0;
			}
			else {
				f.close();
			}
		}
		else { // make new hist-file
			if (!f.fail()) {
				std::cout << "! Wrong: " << temp << " already exists." << std::endl;
				f.close();
				if (argc == 1)
					system("pause");
				return 0;
			}
		}
	}
	// ====== initial file reading finished ======


	std::vector<OBJ> objToDo{};
	for (strType& molName : MolTodo) { // import data from mol file
		strType molPATH = ToDoFolderPath + molName + ".mol";
		OBJ obj0{};
		std::vector<OBJ> objList{};
		int intMsg = obj0.importData(molPATH, objList);
		if (intMsg != 1) {
			std::cout << "! Wrong: mol-file error." << std::endl;
			if (argc == 1)
				system("pause");
			return 0;
		}
		objToDo.push_back(objList[0]);
	}
	for (size_t j = 0; j < objToDo.size(); ++j) {
		// ====== go through step1 ======
		long long int Nstep1 = Nstep1vec[j];
		if (Nstep1 == 0) {
			continue;
		}
		else {
			OBJ obj = objToDo[j];
			unsigned long long int Mt = obj.calTotalM(); // total number of possible cutting schemes
			if (Nstep1 == -1) {
				Nstep1 = static_cast<long long int> (Mt / 100);
				if (Nstep1 < Nstep1_Min) {
					Nstep1 = Nstep1_Min;
				}
			}

			std::chrono::time_point<std::chrono::steady_clock> startT, step1T, endT;
			startT = std::chrono::steady_clock::now();
			std::chrono::duration<double> running_time;
			double ElapsedTime;
			std::cout << std::endl << "Please wait..." << std::endl << std::endl;
			std::cout << " ." << MolTodo[j] << ":" << std::endl;
			std::cout << "  " << Nstep1 << " steps to try in step 1..." << std::endl;
			step1T = std::chrono::steady_clock::now();
			// Monte Carlo step 1: generate fragments histogram
			std::vector<std::tuple<strType, strType, unsigned long long int>> histItem{};
			FRAGHIST fragHist{};
			histItem = obj.histItemOrigMol();
			fragHist.addOnce(histItem);
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
			std::cout << std::endl << std::endl;
			fragHist.getHistOverall(Mt);
			endT = std::chrono::steady_clock::now();
			running_time = endT - startT;
			ElapsedTime = running_time.count();

			// output historam
			strType outPATH_histWhole{};
			outPATH_histWhole = ToDoFolderPath + MolTodo[j] + "_histWhole.txt";
			outputFile_histWhole(obj, fragHist, outPATH_histWhole,
				Nstep1, ElapsedTime);
		}
	}
	std::cout << std::endl
		<< "===============================================" << std::endl
		<< "====== Step 1 is done for all molecules. ======" << std::endl
		<< "===============================================" << std::endl;


	// ==================================================================
	// ========== step 2: from here combine existing histograms =========
	std::cout << std::endl
		<< "====== Step 2 started. ======" << std::endl;

	strType molPATH = ToDoFolderPath + MolTodo[0] + ".mol";
	OBJ obj{};
	std::vector<OBJ> objList{};
	obj.importData(molPATH, objList);
	obj = objList[0];
	FRAGHIST_EXTERNAL histWhole{};
	strType fileName = ToDoFolderPath + MolTodo[0] + "_histWhole.txt";
	int count = histWhole.getHistFromFile(fileName, 0);
	std::cout << "File 1 (" << MolTodo[0] <<") InChI counts: " << count << std::endl;

	for (int i = 1; i < MolTodo.size(); ++i) {
		molPATH = ToDoFolderPath + MolTodo[i] + ".mol";
		int bondIDmax1 = obj.addMol(molPATH);

		fileName = ToDoFolderPath + MolTodo[i] + "_histWhole.txt";
		count = histWhole.getHistFromFile(fileName, bondIDmax1);
		std::cout << "File " << i+1 << " (" << MolTodo[i] << ") InChI counts: " << count << std::endl;
	}
	FRAGHIST hist0 = histWhole.convertToFRAGHIST();
	// histogram info obtained
	// ====== step 2 calculation ======
	strType PATH_AllPaths = ToDoFolderPath + "Tree_AllPaths.txt";
	int Nmol = static_cast<int>(MolTodo.size());
	// ============ mian ============
	if (Nstep2 == -1) {
		Nstep2 = Nstep2_default;
	}
	std::vector<PATHWAY_ASSEMBLY> PAvec = obj.collectPaths(hist0, Nmol, Nstep2);
	outputAllPaths(PAvec, PATH_AllPaths);
	outputTree(PAvec, ToDoFolderPath + "Tree_byRepeats.txt", -1, "nRepeats");
	outputTree(PAvec, ToDoFolderPath + "Tree_bySize.txt", -1, "size");
	std::cout << "Pathway distribution:" << std::endl;
	for (size_t i = 1; i < PAvec.size(); ++i) {
		if (!PAvec[i].paths.empty()) {
			std::cout << i << ":" << PAvec[i].paths.size() << std::endl;
		}
	}
	// ============ mian ============

	
	if (argc == 1)
		system("pause");
	return 1;
}
