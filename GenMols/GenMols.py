# To generated molecules from an assembly pool

import random
import numpy as np
import math
import Reassembler as ra # import our Reassembler
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def pickTwo(ToBeComb):
    """Pick two fragments randomly. ToBeComb records fragments indices."""
    picked, remain = [], []
    if len(ToBeComb) <= 1:
        return picked, remain
    else:
        temp = np.random.choice(len(ToBeComb), 2, replace=False)
        picked = [ToBeComb[temp[0]], ToBeComb[temp[1]]]
        for i in range(len(ToBeComb)):
            if i != temp[0] and i != temp[1]:
                remain.append(ToBeComb[i]) # get what index remained
        return picked, remain


def getNumAtom(formu, ele):
    """Calculate the number of elements (ele) in the molecule (represented by its formula)"""
    idx1 = formu.find(ele)
    if idx1 == -1: # no ele found
        return 0
    else:
        temp = formu[ idx1+len(ele) :]
        if len(temp) == 0: # e.g., C20H25NO, here #O = 1
            return 1
        else:
            try:
                idx2 = temp.find(next(filter(str.isalpha, temp)))
            except StopIteration:
                return int(temp) # when temp is the last digit
            nstr = temp[ : idx2]
            if len(nstr) == 0: # e.g., C20H25NO, here #N = 1
                return 1
            else:
                return int(nstr)


def DegUnsat(mol):
    """Calculate Degree of Unsaturation of mol"""
    formu = rdMolDescriptors.CalcMolFormula(mol)
    nC = getNumAtom(formu, 'C')
    nN = getNumAtom(formu, 'N')
    nX = getNumAtom(formu, 'F') + getNumAtom(formu, 'Cl')          + getNumAtom(formu, 'Br') + getNumAtom(formu, 'I')
    nH = getNumAtom(formu, 'H')
    dou = (2*nC + 2 + nN - nX - nH) / 2 # a well-defined formula
    return dou


def GenerateMols(
    NmolNeeded,
    PoolFile,
    OutputInchiFile,
    OutputFigPath,
    mwMin,
    mwMax,
    DoUMin,
    DoUMax,
    oneAtomWeight = 12,
    mwDelta = 0.1
):
    """The main function that generates new molecules from an assembly pool.
    Args:
        NmolNeeded (int): how many molecules needed to be generated.
        PoolFile (str): path to assembly pool inchis.
        OutputInchiFile (str): file name where the output inchis will be written into.
        OutputFigPath (str): path of a folder where the pictures of the newly-generated molecules will be put into.
        mwMin (int): minimum molecular weight of newly-generated molecules.
        mwMax (int): maximum molecular weight of newly-generated molecules.
        DoUMin (int): minimum Degree of Unsaturation of newly-generated molecules.
        DoUMax (int): maximum Degree of Unsaturation of newly-generated molecules.
        oneAtomWeight (int): an approximate molecular weight lost when an atom is thrown away when two fragments are combined.
        mwDelta (float): give a range of the molecular weight that could be relaxed.
    """
    frag = []
    with open(PoolFile) as f:
        for inc in f:
            if inc[:5] == 'InChI':
                inc = inc.replace('\n', '')
                frag.append(Chem.MolFromInchi(inc))
    Nfrag = len(frag)
    print('# fragments in the assembly pool:', Nfrag, flush=True)
    
    mwMinTry = mwMin * (1 - mwDelta)
    mwMaxTry = mwMax * (1 + mwDelta)
    n = 0
    newMols = []
    nFiltered = 0

    while n < NmolNeeded:
        try:
            mw = oneAtomWeight
            idxList = []
            idx = []
            while True:
            # to obtain idxList that contains a list of (list of fragments)
            # each (list of fragments) is within the molecular weight range [mwMinTry, mwMaxTry]
            # so later, we randomly choose one (list of fragments) is a proper set of fragments that when they're combined, the weight is within the requred range.
                idxThis = random.randint(0, Nfrag-1)
                mw += rdMolDescriptors.CalcExactMolWt(frag[idxThis], onlyHeavy=True)
                if mw < mwMinTry:
                    idx.append(idxThis)
                    mw -= oneAtomWeight
                else:
                    if mw <= mwMaxTry:
                        idx.append(idxThis)
                        idxList.append(idx.copy())
                        mw -= oneAtomWeight
                    else:
                        break

            ContinueSearch = False
            if len(idxList) == 0:
                print('warning: mwMin and mwMax may be too close.', flush=True)
            else:
                ToBeComb = idxList[random.randint(0, len(idxList)-1)]
                # randomly choose a (list of fragments)

                ifrag = Nfrag
                while len(ToBeComb) > 1:
                # combine fragments until only a big one exists.
                # randomly pick two fragments from ToBeComb and combine them, and then put it back into ToBeComb; Then randomly pick two from ToBeComb and repeat.
                    picked, remain = pickTwo(ToBeComb)
                    M1 = frag[picked[0]]
                    M2 = frag[picked[1]]

                    if M1.GetNumBonds() == 1 or M2.GetNumBonds() == 1:
                        # when either has only one bond, one one atom can be joined ("overlapped")
                        newmol = ra.assemble(M1, M2, 1)
                    else:
                        # otherwise, decide whether overlap 1 or 2 atoms
                        newmol = ra.assemble(M1, M2, random.randint(1, 2))

                    if newmol is None:
                        # it could happen that combination fails
                        ContinueSearch = True
                        break
                    else:
                        if ifrag >= len(frag):
                            frag.append(newmol)
                        else:
                            frag[ifrag] = newmol
                        remain.append(ifrag)
                        ifrag += 1
                        ToBeComb = remain
                if ContinueSearch:
                    continue
                newMolecule = frag[ToBeComb[0]]

                # consider Degree of Unsaturation (DoU)
                dou = DegUnsat(newmol)
                if dou > DoUMax:
                    print('>', end='', flush=True)
                    continue
                elif dou < DoUMin:
                    print('<', end='', flush=True)
                    # if DoU is smaller than DoUMin, then use operation Origami() to make rings (details in Reassembler.py)
                    # one Origami() operation, DoU is increased by 1.
                    nOrigami = random.randint(DoUMin-int(dou), int(DoUMax-dou))
                    # randomly decide how many Origami() operation will be done
                    try:
                        for i in range(nOrigami):
                            newMolecule = random.choice(ra.origami(newMolecule))
                    except:
                        print('-', end='', flush=True)
                        continue
                else:
                    pass


                # apply filterMol() which is commonly used to filter obviously impossible molecules
                newMolecule = ra.filterMol(newMolecule)
                if newMolecule is None:
                    print(':', end='', flush=True)
                    nFiltered += 1
                    continue

                # apply confilter() which filters molecules that do not have valid conformations
                newMolecule = ra.confilter(newMolecule)
                if newMolecule is None:
                    print(':', end='', flush=True)
                    nFiltered += 1
                    continue

                # then, a new molecule is found
                n += 1
                print(n, end=', ', flush=True)
                newMols.append(Chem.MolToInchi(newMolecule))
        except:
            continue


    f = open(OutputInchiFile, 'w')
    i = 0
    for inc in newMols:
        i += 1
        f.write(inc + '\n')
        f.write('--- Molecule ' + str(i) + '\n')
    f.close()
    print('\n', 'nFiltered =', nFiltered, flush=True)

    # visualize, generate figures
    ra.printer(OutputInchiFile, OutputFigPath)


def run(NmolNeeded, PoolFile = 'Pool.txt', OutputInchiFile = 'newMols.txt', OutputFigPath = 'MolsFig', mwMin = 281, mwMax = 368, DoUMin = 9, DoUMax = 12, oneAtomWeight = 12, mwDelta = 0.1):

#    NmolNeeded = 1000
#    PoolFile = 'Pool.txt'
#    OutputFile = 'newMols.txt'
#    OutputFigPath = 'MolsFig'
#    mwMin = 281 # minimum molecular weight of the 6 natural opiates we used
#    mwMax = 368 # maximum molecular weight of the 6 natural opiates we used
#    DoUMin = 9 # minimum Degree of Unsaturation of the 6 natural opiates we used
#    DoUMax = 12 # maximum Degree of Unsaturation of the 6 natural opiates we used
#    oneAtomWeight = 12 # choose for Carbon atom
#    mwDelta = 0.1

    GenerateMols(NmolNeeded, PoolFile, OutputInchiFile, OutputFigPath, mwMin, mwMax, DoUMin, DoUMax, oneAtomWeight, mwDelta)
    print('Done.')
