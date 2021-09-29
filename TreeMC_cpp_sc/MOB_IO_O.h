#pragma once
#include "MOB_Def.h"
#include "MOB.h"


int outputFile(MOL_BOND& mol, PATHWAY_ASSEMBLY& PA, int stage,
	strType& outPATH, strType& outPATH_hist,
	long long int Nstep1, long long int Nstep2,
	long long int nTotalTry, double ElapsedTime);


//void outputFile_hist(MOL_BOND& mol, FRAGHIST& fragHist,
//	strType& outPATH_hist,
//	long long int Nstep1, double ElapsedTime);
void outputFile_histWhole(MOL_BOND& mol, FRAGHIST& fragHist,
	strType& outPATH_histWhole,
	long long int Nstep1, double ElapsedTime);


void outputAllPaths(std::vector<PATHWAY_ASSEMBLY>& PAvec, strType FilePath);
void outputTree(std::vector<PATHWAY_ASSEMBLY>& PAvec,
	strType FilePath, int ThresholdPA, strType orderBy);


//void outputConsole(MOL_BOND& mol, PATHWAY_ASSEMBLY& PA,
//	strType& outPATH, strType& outPATH_hist,
//	long long int Nstep1, long long int Nstep2,
//	long long int nTotalTry, double ElapsedTime);


//void printInstruction();
