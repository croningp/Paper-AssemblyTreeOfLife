#include "MOB_Inchi.h"
#include "inchi_api.h"
#include <unordered_map>
// Note that I use #include "MOB_Inchi_Appendix.h" in the text

/* =========================================== */
strType getInchi(MOL_BOND& mol) {
	int i = -1, k;
	AT_NUM nAtoms = static_cast<AT_NUM>(mol.Atoms.size());
	// 1st part of the function is to otbain the inchi_Atom class that 
	// InChI API needs for the input of InChI string generator.
	inchi_Atom* inchiMol = new inchi_Atom[nAtoms]; // defined in external  InChI API
	char AtomicSym[ATOM_EL_LEN] = { '\0','\0','\0','\0','\0','\0' };
	std::unordered_map<int, int> tempID;

	for (auto& atom : mol.Atoms) {
		i++;
		#include "MOB_Inchi_Appendix.h" // Note here
		strcpy_s(inchiMol[i].elname, AtomicSym);
		inchiMol[i].x = 0; // no coordinates
		inchiMol[i].y = 0; // no coordinates
		inchiMol[i].z = 0; // no coordinates
		inchiMol[i].num_iso_H[0] = -1; // implicit hydrogens based on valence
		inchiMol[i].num_iso_H[1] = 0; // no isotope
		inchiMol[i].num_iso_H[2] = 0; // no isotope
		inchiMol[i].num_iso_H[3] = 0; // no isotope
		inchiMol[i].isotopic_mass = 0; // no isotopic mass
		inchiMol[i].radical = 0; // no radical
		inchiMol[i].charge = 0; // zero charge

		k = -1;
		for (int& j : atom.second.BondsLinked) {
			k++;
			if (atom.first == mol.Bonds[j].Atom1id) {
				inchiMol[i].neighbor[k] = static_cast<AT_NUM>(mol.Bonds[j].Atom2id);
			}
			else { inchiMol[i].neighbor[k] = static_cast<AT_NUM>(mol.Bonds[j].Atom1id); }
			inchiMol[i].bond_type[k] = mol.Bonds[j].BondType;
			inchiMol[i].bond_stereo[k] = 3; // 3 means unknown
		}
		inchiMol[i].num_bonds = k + 1;
		// tempID records the atom's orignal index in MOL, 
		// associated with the temporary index here
		tempID[atom.first] = i;
	}
	for (int i = 0; i < nAtoms; ++i) {
		for (int j = 0; j < inchiMol[i].num_bonds; ++j) {
			inchiMol[i].neighbor[j] = tempID[inchiMol[i].neighbor[j]];
		}
	}


	// inchi_Input class is the data type the InChI generator needs.
	inchi_Input inchi_input_data; // defined in external  InChI API
	inchi_input_data.atom = inchiMol;
	inchi_input_data.num_atoms = nAtoms;
	inchi_input_data.num_stereo0D = 0; // no stereo information
	inchi_input_data.stereo0D = NULL;
	// "fixed H" option is on, this is to make sure the generated InChI 
	// distinguishes tautomerism.
	char opt[] = "/FixedH";
	inchi_input_data.szOptions = opt; // no options
		//char opt = ' ';
		//inchi_input_data.szOptions = &opt; // no options

	inchi_Output inchi_output_data;
	GetINCHI(&inchi_input_data, &inchi_output_data); // get the inchi (the "human readable" string)


	strType output_inchi_string = inchi_output_data.szInChI;
	//std::string inchi_log = inchi_output_data.szLog;

	// This is to calculate InChIKey, but we don't need it here.
	//char* szINCHIKey = new char[28];
	//char* szExtra1 = new char[65];
	//char* szExtra2 = new char[65];
	//GetINCHIKeyFromINCHI(inchi_output_data.szInChI, 1, 1, szINCHIKey, szExtra1, szExtra2); // get the inchi key (hashed fixed_length version of the inchi string)
	//std::string output_key(szINCHIKey);

	FreeStdINCHI(&inchi_output_data); // clean up pointers etc
	delete[] inchiMol; // delete the fragment_atoms object
	//delete[] szINCHIKey; delete[] szExtra1; delete[] szExtra2;

	//if (key) return output_key;
	//else return output_inchi_string;

	return output_inchi_string;
}