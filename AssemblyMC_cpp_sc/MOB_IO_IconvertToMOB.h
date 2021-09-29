/* It defines functions that import data (IO_I_readMolFile::MolData)
   to MOL that is used in the acutal calculation.
*/
#pragma once
#include "MOB_Def.h"
#include "MOB_IO_IreadMolFile.h"
#include "MOB.h"


// It transforms molfile::MolData to MOL
// that we use in the actual calculations.
MOL_BOND getMOL_from(MOB_IO_IreadMolFile::MolData& moldata, int nfrag0);


// Associate Atomic Number with Atomic Symbol, and vice versa.
int getAtomicNum(strType& AtomSym);
strType getAtomicSym2Char(int AtomicNum);
