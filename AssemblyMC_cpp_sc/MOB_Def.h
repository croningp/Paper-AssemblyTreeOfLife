/* It declares definitions and constans that all functions need. */
#pragma once
#include <vector>
#include <string>
typedef std::vector<int> intVec;
typedef std::vector< std::vector<int> > intVec2D;
typedef std::string strType;


const std::string versionStr = "Monte Carlo Calculating Pathway Assembly (in bonds). v1.0";


const std::string Err_Msg_900 = "! Mol-file read successfully, but it contains no valid molecules.";
const std::string Err_Msg_902 = "! Mol-file read successfully, but it contains >= 2 valid molecules. This program only deals with single molecule.";


// Ntry_freeze: how many times GunShot generates no-change, then stop
// if Ntry_freeze = -1, then check bond.size() * Ntry_freeze_amplifier times
const int Ntry_freeze_amplifier = 2;


const int Nstep1_Min = 100000; //100,000
const int Nstep2_SaveWindow = 10000; //10,000, must be 10*
