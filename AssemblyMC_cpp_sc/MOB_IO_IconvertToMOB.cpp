#include "MOB_IO_IconvertToMOB.h"


/* =========================================== */
int getAtomicNum(strType& AtomSym) {
	int i = 0;
	if (AtomSym == "H") { i = 1; return i; }
	if (AtomSym == "He") { i = 2; return i; }
	if (AtomSym == "Li") { i = 3; return i; }
	if (AtomSym == "Be") { i = 4; return i; }
	if (AtomSym == "B") { i = 5; return i; }
	if (AtomSym == "C") { i = 6; return i; }
	if (AtomSym == "N") { i = 7; return i; }
	if (AtomSym == "O") { i = 8; return i; }
	if (AtomSym == "F") { i = 9; return i; }
	if (AtomSym == "Ne") { i = 10; return i; }

	if (AtomSym == "Na") { i = 11; return i; }
	if (AtomSym == "Mg") { i = 12; return i; }
	if (AtomSym == "Al") { i = 13; return i; }
	if (AtomSym == "Si") { i = 14; return i; }
	if (AtomSym == "P") { i = 15; return i; }
	if (AtomSym == "S") { i = 16; return i; }
	if (AtomSym == "Cl") { i = 17; return i; }
	if (AtomSym == "Ar") { i = 18; return i; }

	if (AtomSym == "K") { i = 19; return i; }
	if (AtomSym == "Ca") { i = 20; return i; }
	if (AtomSym == "Sc") { i = 21; return i; }
	if (AtomSym == "Ti") { i = 22; return i; }
	if (AtomSym == "V") { i = 23; return i; }
	if (AtomSym == "Cr") { i = 24; return i; }
	if (AtomSym == "Mn") { i = 25; return i; }
	if (AtomSym == "Fe") { i = 26; return i; }
	if (AtomSym == "Co") { i = 27; return i; }
	if (AtomSym == "Ni") { i = 28; return i; }
	if (AtomSym == "Cu") { i = 29; return i; }
	if (AtomSym == "Zn") { i = 30; return i; }
	if (AtomSym == "Ga") { i = 31; return i; }
	if (AtomSym == "Ge") { i = 32; return i; }
	if (AtomSym == "As") { i = 33; return i; }
	if (AtomSym == "Se") { i = 34; return i; }
	if (AtomSym == "Br") { i = 35; return i; }
	if (AtomSym == "Kr") { i = 36; return i; }

	if (AtomSym == "Rb") { i = 37; return i; }
	if (AtomSym == "Sr") { i = 38; return i; }
	if (AtomSym == "Y") { i = 39; return i; }
	if (AtomSym == "Zr") { i = 40; return i; }
	if (AtomSym == "Nb") { i = 41; return i; }
	if (AtomSym == "Mo") { i = 42; return i; }
	if (AtomSym == "Tc") { i = 43; return i; }
	if (AtomSym == "Ru") { i = 44; return i; }
	if (AtomSym == "Rh") { i = 45; return i; }
	if (AtomSym == "Pd") { i = 46; return i; }
	if (AtomSym == "Ag") { i = 47; return i; }
	if (AtomSym == "Cd") { i = 48; return i; }
	if (AtomSym == "In") { i = 49; return i; }
	if (AtomSym == "Sn") { i = 50; return i; }
	if (AtomSym == "Sb") { i = 51; return i; }
	if (AtomSym == "Te") { i = 52; return i; }
	if (AtomSym == "I") { i = 53; return i; }
	if (AtomSym == "Xe") { i = 54; return i; }

	if (AtomSym == "Cs") { i = 55; return i; }
	if (AtomSym == "Ba") { i = 56; return i; }
	if (AtomSym == "La") { i = 57; return i; }

	// transition metal
	if (AtomSym == "Ce") { i = 58; return i; }
	if (AtomSym == "Pr") { i = 59; return i; }
	if (AtomSym == "Nd") { i = 60; return i; }
	if (AtomSym == "Pm") { i = 61; return i; }
	if (AtomSym == "Sm") { i = 62; return i; }
	if (AtomSym == "Eu") { i = 63; return i; }
	if (AtomSym == "Gd") { i = 64; return i; }
	if (AtomSym == "Tb") { i = 65; return i; }
	if (AtomSym == "Dy") { i = 66; return i; }
	if (AtomSym == "Ho") { i = 67; return i; }
	if (AtomSym == "Er") { i = 68; return i; }
	if (AtomSym == "Tm") { i = 69; return i; }
	if (AtomSym == "Yb") { i = 70; return i; }
	if (AtomSym == "Lu") { i = 71; return i; }

	if (AtomSym == "Hf") { i = 72; return i; }
	if (AtomSym == "Ta") { i = 73; return i; }
	if (AtomSym == "W") { i = 74; return i; }
	if (AtomSym == "Re") { i = 75; return i; }
	if (AtomSym == "Os") { i = 76; return i; }
	if (AtomSym == "Ir") { i = 77; return i; }
	if (AtomSym == "Pt") { i = 78; return i; }
	if (AtomSym == "Au") { i = 79; return i; }
	if (AtomSym == "Hg") { i = 80; return i; }
	if (AtomSym == "Tl") { i = 81; return i; }
	if (AtomSym == "Pb") { i = 82; return i; }
	if (AtomSym == "Bi") { i = 83; return i; }
	if (AtomSym == "Po") { i = 84; return i; }
	if (AtomSym == "At") { i = 85; return i; }
	if (AtomSym == "Rn") { i = 86; return i; }

	if (AtomSym == "Fr") { i = 87; return i; }
	if (AtomSym == "Ra") { i = 88; return i; }
	if (AtomSym == "Ac") { i = 89; return i; }

	// transition metal
	if (AtomSym == "Th") { i = 90; return i; }
	if (AtomSym == "Pa") { i = 91; return i; }
	if (AtomSym == "U") { i = 92; return i; }
	if (AtomSym == "Np") { i = 93; return i; }
	if (AtomSym == "Pu") { i = 94; return i; }
	if (AtomSym == "Am") { i = 95; return i; }
	if (AtomSym == "Cm") { i = 96; return i; }
	if (AtomSym == "Bk") { i = 97; return i; }
	if (AtomSym == "Cf") { i = 98; return i; }
	if (AtomSym == "Es") { i = 99; return i; }
	if (AtomSym == "Fm") { i = 100; return i; }
	if (AtomSym == "Md") { i = 101; return i; }
	if (AtomSym == "No") { i = 102; return i; }
	if (AtomSym == "Lr") { i = 103; return i; }

	if (AtomSym == "Rf") { i = 104; return i; }
	if (AtomSym == "Ha") { i = 105; return i; }
	if (AtomSym == "Sg") { i = 106; return i; }

	return i;
}


/* =========================================== */
// All atomic symbols have size 2 even "H ",
// for convenience of output.
strType getAtomicSym2Char(int AtomicNum) {
	switch (AtomicNum) {
	case 1: return "H ";
	case 2: return "He";
	case 3: return "Li";
	case 4: return "Be";
	case 5: return "B ";
	case 6: return "C ";
	case 7: return "N ";
	case 8: return "O ";
	case 9: return "F ";
	case 10: return "Ne";

	case 11: return "Na";
	case 12: return "Mg";
	case 13: return "Al";
	case 14: return "Si";
	case 15: return "P ";
	case 16: return "S ";
	case 17: return "Cl";
	case 18: return "Ar";

	case 19: return "K ";
	case 20: return "Ca";
	case 21: return "Sc";
	case 22: return "Ti";
	case 23: return "V ";
	case 24: return "Cr";
	case 25: return "Mn";
	case 26: return "Fe";
	case 27: return "Co";
	case 28: return "Ni";
	case 29: return "Cu";
	case 30: return "Zn";
	case 31: return "Ga";
	case 32: return "Ge";
	case 33: return "As";
	case 34: return "Se";
	case 35: return "Br";
	case 36: return "Kr";

	case 37: return "Rb";
	case 38: return "Sr";
	case 39: return "Y ";
	case 40: return "Zr";
	case 41: return "Nb";
	case 42: return "Mo";
	case 43: return "Tc";
	case 44: return "Ru";
	case 45: return "Rh";
	case 46: return "Pd";
	case 47: return "Ag";
	case 48: return "Cd";
	case 49: return "In";
	case 50: return "Sn";
	case 51: return "Sb";
	case 52: return "Te";
	case 53: return "I ";
	case 54: return "Xe";

	case 55: return "Cs";
	case 56: return "Ba";
	case 57: return "La";

		// transition metal
	case 58: return "Ce";
	case 59: return "Pr";
	case 60: return "Nd";
	case 61: return "Pm";
	case 62: return "Sm";
	case 63: return "Eu";
	case 64: return "Gd";
	case 65: return "Tb";
	case 66: return "Dy";
	case 67: return "Ho";
	case 68: return "Er";
	case 69: return "Tm";
	case 70: return "Yb";
	case 71: return "Lu";

	case 72: return "Hf";
	case 73: return "Ta";
	case 74: return "W ";
	case 75: return "Re";
	case 76: return "Os";
	case 77: return "Ir";
	case 78: return "Pt";
	case 79: return "Au";
	case 80: return "Hg";
	case 81: return "Tl";
	case 82: return "Pb";
	case 83: return "Bi";
	case 84: return "Po";
	case 85: return "At";
	case 86: return "Rn";

	case 87: return "Fr";
	case 88: return "Ra";
	case 89: return "Ac";

		// transition metal
	case 90: return "Th";
	case 91: return "Pa";
	case 92: return "U ";
	case 93: return "Np";
	case 94: return "Pu";
	case 95: return "Am";
	case 96: return "Cm";
	case 97: return "Bk";
	case 98: return "Cf";
	case 99: return "Es";
	case 100: return "Fm";
	case 101: return "Md";
	case 102: return "No";
	case 103: return "Lr";

	case 104: return "Rf";
	case 105: return "Ha";
	case 106: return "Sg";

	default: return "? ";
	}
}


/* =========================================== */
// It transforms molfile::MolData to MOL 
// that is used in the acutal calculation.
MOL_BOND getMOL_from(MOB_IO_IreadMolFile::MolData& moldata, int nfrag0) {
	intVec AtomsList{};
	for (strType& atomSym : moldata.atoms) {
		AtomsList.push_back(getAtomicNum(atomSym));
	}
	MOL_BOND mol(AtomsList, moldata.bonds, nfrag0);

	return mol;
}
