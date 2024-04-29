#!/usr/bin/python3
""" Python script to generate all atomic density matrices for basis sets in the NWChem library 
    with NWChem calculations

# Usage: to run this script, 1) If you only want to generate atomic density matrice for 
                                several basis sets, you can comment out the line 
                                "for basis_name in basis_elements:" and uncomment the line
                                "for basis_name in small_basis_elements:", then copy the
                                corresponding key-value pairs in the basis_elements dictionary
                                to the small_basis_elements dictionary; 
                             1') Or you just want to generate atomic density matrices for all
                                available basis sets, you can comment out the line 
                                "for basis_name in small_basis_elements:" and uncomment the line
                                "for basis_name in basis_elements:";
#                            2) Modify the line "subprocess.run(["./mynw.sh", 'tmp'])" below,
#                               make it adapted to the path to your NWChem script or executable. 
#                               Note here the input file is "tmp.nw" and the output file is set 
#                               to "tmp.nwo" by default.
#                            3) Make two directories "scratch" and "perm" under the current
#                               directory to hold the scratch and permanent files of NWChem
#                               calculations.
#                            4) Simply run the script in the command line as: ./batch_nwchem_dm.py
#                               The script will try to map the basis set in your choice of 
#                               dictionary to this list of NWChemEx basis sets. If the mapping is
#                               successful, atomic density matrix data files with the basis set
#                               names in the NWChemEx convention will be generated.
#
# The script write a NWChem input file (tmp.nw) for each atom in the value list corresponding to
# the basis name key, and call NWChem to run a simple SAD guess generation calculaton (noscf).
#
# After the calculation the SAD guess (atomic density matrix) is extracted and reformmated, and
# stored accumulatively in a "basis_name.dat" file (e. g., STO-3G.dat).
#
# The atomic density matrices are stored as text blocks representing the 1-D flatten version of
# the corresponding 2-D matrices.
# This script works only for basis sets defined in the NWChem library (see the dictionary 
# basis_elements below). For the same calculation with user defined basis sets, please use
# another script genbase_nwchem_dm.py.
"""

import numpy as np
import subprocess
import os


def atom_sym2num(sym):
    elements = {
        'H': 1,
        'He': 2,
        'Li': 3,
        'Be': 4,
        'B': 5,
        'C': 6,
        'N': 7,
        'O': 8,
        'F': 9,
        'Ne': 10,
        'Na': 11,
        'Mg': 12,
        'Al': 13,
        'Si': 14,
        'P': 15,
        'S': 16,
        'Cl': 17,
        'Ar': 18,
        'K': 19,
        'Ca': 20,
        'Sc': 21,
        'Ti': 22,
        'V': 23,
        'Cr': 24,
        'Mn': 25,
        'Fe': 26,
        'Co': 27,
        'Ni': 28,
        'Cu': 29,
        'Zn': 30,
        'Ga': 31,
        'Ge': 32,
        'As': 33,
        'Se': 34,
        'Br': 35,
        'Kr': 36,
        'Rb': 37,
        'Sr': 38,
        'Y': 39,
        'Zr': 40,
        'Nb': 41,
        'Mo': 42,
        'Tc': 43,
        'Ru': 44,
        'Rh': 45,
        'Pd': 46,
        'Ag': 47,
        'Cd': 48,
        'In': 49,
        'Sn': 50,
        'Sb': 51,
        'Te': 52,
        'I': 53,
        'Xe': 54,
        'Cs': 55,
        'Ba': 56,
        'La': 57,
        'Ce': 58,
        'Pr': 59,
        'Nd': 60,
        'Pm': 61,
        'Sm': 62,
        'Eu': 63,
        'Gd': 64,
        'Tb': 65,
        'Dy': 66,
        'Ho': 67,
        'Er': 68,
        'Tm': 69,
        'Yb': 70,
        'Lu': 71,
        'Hf': 72,
        'Ta': 73,
        'W': 74,
        'Re': 75,
        'Os': 76,
        'Ir': 77,
        'Pt': 78,
        'Au': 79,
        'Hg': 80,
        'Tl': 81,
        'Pb': 82,
        'Bi': 83,
        'Po': 84,
        'At': 85,
        'Rn': 86,
        'Fr': 87,
        'Ra': 88,
        'Ac': 89,
        'Th': 90,
        'Pa': 91,
        'U': 92,
        'Np': 93,
        'Pu': 94,
        'Am': 95,
        'Cm': 96,
        'Bk': 97,
        'Cf': 98,
        'Es': 99,
        'Fm': 100,
        'Md': 101,
        'No': 102,
        'Lr': 103,
        'Rf': 104,
        'Db': 105,
        'Sg': 106,
        'Bh': 107,
        'Hs': 108,
        'Mt': 109,
        'Ds': 110,
        'Rg': 111,
        'Cn': 112,
        'Nh': 113,
        'Fl': 114,
        'Mc': 115,
        'Lv': 116,
        'Ts': 117,
        'Og': 118
    }

    try:
        return elements[sym]
    except:
        print("Cannot find the atomic number for", sym, "!")


# To get the atomic symbol from the corresponding atomic number
#print(list(elements.keys())[list(elements.values()).index(3)])

# The following dictionary contains the information of the basis set names and
# the atoms for which the corresponding basis sets are available
# Info refomated from https://nwchemgit.github.io/AvailableBasisSets.html
# Basis set names are in the convention of NWChem

basis_elements = {
    '3-21++G': '  H Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca ',
    '3-21++G*': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    '3-21G':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs ',
    '3-21G*': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    '3-21GSP': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    '4-22GSP': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    '4-31G': '  H He Li Be B C N O F Ne P S Cl ',
    '6-31++G': '  H Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca ',
    '6-31++G*': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca ',
    '6-31++G**': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca ',
    '6-31+G*': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca ',
    '6-311++G(2d,2p)': '  H Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca ',
    '6-311++G(3df,3pd)': '  H Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    '6-311++G**': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca ',
    '6-311+G*': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca ',
    '6-311G':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Ga Ge As Se Br Kr I ',
    '6-311G(2df,2pd)': '  H He Li Be B C N O F Ne K Ca ',
    '6-311G*':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Ga Ge As Se Br Kr I ',
    '6-311G**':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Ga Ge As Se Br Kr I ',
    '6-31G':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    '6-31G-Blaudeau': '  K Ca ',
    '6-31G(3df,3pd)': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    '6-31G*':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    '6-31G*-Blaudeau': '  K Ca ',
    '6-31G**':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'Ahlrichs pVDZ':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'Ahlrichs TZV':
    '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'Ahlrichs_VDZ':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'Ahlrichs VTZ':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'ANO-RCC':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm ',
    'apr-cc-pV(Q+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'aug-cc-pCV5Z': '  B C N O F Ne Al Si P S Cl Ar ',
    'aug-cc-pCVDZ': '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'aug-cc-pCVQZ': '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'aug-cc-pCV(T+d)Z': '  Al Si P S Cl ',
    'aug-cc-pCVTZ': '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'aug-cc-pV(5+d)Z': '  Al Si P S Cl Ar ',
    'aug-cc-pV5Z':
    '  H He B C N O F Ne Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-cc-pV(6+d)Z': '  Al Si P S Cl Ar ',
    'aug-cc-pV6Z': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-cc-pV(D+d)Z': '  Al Si P S Cl Ar ',
    'aug-cc-pVDZ':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-cc-pV(Q+d)Z': '  Al Si P S Cl Ar ',
    'aug-cc-pVQZ':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-cc-pV(T+d)Z': '  Al Si P S Cl Ar ',
    'aug-cc-pVTZ':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-cc-pVTZ-J':
    '  H B C N O F Al Si P S Cl Sc Ti V Cr Mn Fe Co Ni Cu Zn Se ',
    'aug-cc-pwCV5Z': '  B C N O F Ne Al Si P S Cl Ar ',
    'aug-cc-pwCV5Z-NR': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'aug-cc-pwCVDZ': '  B C N O F Ne Al Si P S Cl Ar ',
    'aug-cc-pwCVQZ': '  B C N O F Ne Al Si P S Cl Ar Br ',
    'aug-cc-pwCVQZ-NR': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'aug-cc-pwCVTZ': '  B C N O F Ne Al Si P S Cl Ar ',
    'aug-cc-pwCVTZ-NR': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'aug-mcc-pV5Z': '  H ',
    'aug-mcc-pV6Z': '  H ',
    'aug-mcc-pV7Z': '  H ',
    'aug-mcc-pV8Z': '  H ',
    'aug-mcc-pVQZ': '  H ',
    'aug-mcc-pVTZ': '  H ',
    'aug-pc-0':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Ga Ge As Se Br Kr ',
    'aug-pc-1':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-pc-2':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-pc-3':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-pc-4':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-pcJ-0': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-pcJ-0_2006': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-pcJ-1': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-pcJ-1_2006': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-pcJ-2': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-pcJ-2_2006': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-pcJ-3': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-pcJ-3_2006': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-pcJ-4': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-pcJ-4_2006': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-pcS-0': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'aug-pcS-1': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'aug-pcS-2': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'aug-pcS-3': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'aug-pcS-4': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'aug-pV7Z': '  C N O F S ',
    'B2': '  Zn ',
    'Bauschlicher ANO': '  Sc Ti V Cr Mn Fe Co Ni Cu ',
    'Binning/Curtiss SV': '  Ga Ge As Se Br Kr ',
    'Binning/Curtiss SVP': '  Ga Ge As Se Br Kr ',
    'Binning/Curtiss VTZ': '  Ga Ge As Se Br Kr ',
    'Binning/Curtiss VTZP': '  Ga Ge As Se Br Kr ',
    'cc-pCV5Z': '  B C N O F Ne Na Mg Al Si P S Cl Ar Ca ',
    'cc-pCV6Z': '  B C N O F Ne Al Si P S Cl Ar ',
    'cc-pCV6Z(old)': '  O ',
    'cc-pCVDZ': '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ca ',
    'cc-pCVDZ(old)': '  Al Si P S Cl Ar ',
    'cc-pCVQZ': '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ca ',
    'cc-pCVQZ(old)': '  Al Si P S Cl Ar ',
    'cc-pCVTZ': '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ca ',
    'cc-pCVTZ(old)': '  Al Si P S Cl Ar ',
    'cc-pV(5+d)Z': '  Al Si P S Cl Ar ',
    'cc-pV5Z':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'cc-pV(6+d)Z': '  Al Si P S Cl Ar ',
    'cc-pV6Z': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'cc-pV8Z': '  H Ne ',
    'cc-pV9Z': '  Ne ',
    'cc-pV(D+d)Z': '  Al Si P S Cl Ar ',
    'cc-pVDZ':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'cc-pVDZ(seg-opt)':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pV(Q+d)Z': '  Al Si P S Cl Ar ',
    'cc-pVQZ':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'cc-pVQZ(seg-opt)': '  H He C O ',
    'cc-pV(T+d)Z': '  Al Si P S Cl Ar ',
    'cc-pVTZ':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'cc-pVTZ(seg-opt)':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pwCV5Z': '  B C N O F Ne Al Si P S Cl Ar ',
    'cc-pwCV5Z': '  B C N O F Ne Al Si P S Cl Ar ',
    'cc-pwCV5Z-NR': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'cc-pwCVDZ': '  B C N O F Ne Al Si P S Cl Ar ',
    'cc-pwCVQZ': '  B C N O F Ne Al Si P S Cl Ar Br ',
    'cc-pwCVQZ-NR': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'cc-pwCVTZ': '  B C N O F Ne Al Si P S Cl Ar ',
    'cc-pwCVTZ-NR': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'ccemd-2': '  He H Be Li Ne C N F B O Mg Na Ar P Cl Al Si S ',
    'ccemd-3': '  He H Be Li Ne C N F B O Mg Na Ar P Cl Al Si S ',
    'ccJ-pV5Z': '  H He B C N O F Ne ',
    'ccJ-pVDZ': '  H He B C N O F Ne ',
    'ccJ-pVQZ': '  H He B C N O F Ne ',
    'ccJ-pVTZ': '  H He B C N O F Ne ',
    'coemd-2': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'coemd-3': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'coemd-4': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'coemd-ref': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'CVTZ': '  Li Be Na Mg K Ca ',
    'd-aug-cc-pV5Z': '  H He B C N O F Ne ',
    'd-aug-cc-pV6Z': '  H B C N O ',
    'd-aug-cc-pVDZ': '  H He B C N O F Ne ',
    'd-aug-cc-pVQZ': '  H He B C N O F Ne ',
    'd-aug-cc-pVTZ': '  H He B C N O F Ne ',
    'Def2-QZVP':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Def2-QZVPD':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Def2-SVP':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Def2-SVPD':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Def2-TZVP':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Def2-TZVPD':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'dhf-QZVP':
    '  Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'dhf-SV(P)':
    '  Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'dhf-TZVP':
    '  Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Dunning-Hay Double Rydberg': '  B C N O F Ne Al Si P S Cl Ar ',
    'Dunning-Hay Rydberg': '  B C N O F Ne Al Si P S Cl Ar ',
    'DZ+Double Rydberg': '  B C N O F Ne Al Si P S Cl ',
    'DZ+Rydberg': '  B C N O F Ne Al Si P S Cl ',
    'DZ': '  H Li B C N O F Ne Al Si P S Cl ',
    'DZP+Rydberg': '  B C N O F Ne Al Si P S Cl ',
    'DZP': '  H Li B C N O F Ne Al Si P S Cl ',
    'DZQ': '  Y Zr Nb Mo Tc Ru Rh Pd Ag ',
    'DZVP2':
    '  H He Li Be B C N O F Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'DZVP':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe ',
    'Feller_Misc._CVDZ': '  K ',
    'Feller': '  K ',
    'Feller Misc. CVTZ': '  K ',
    'GAMESS PVTZ': '  H Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'GAMESS VTZ': '  H Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'IGLO-II': '  H B C N O F Al Si P S Cl ',
    'IGLO-III': '  H B C N O F Al Si P S Cl ',
    'jul-cc-pV(D+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'jul-cc-pV(Q+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'jul-cc-pV(T+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'jun-cc-pV(D+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'jun-cc-pV(Q+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'jun-cc-pV(T+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'LANL08':
    '  Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi ',
    'LANL08+': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'LANL08d': '  Si P S Cl Ge As Se Br Sn Sb Te I Pb Bi ',
    'LANL08(f)':
    '  Sc Ti V Cr Mn Fe Co Ni Cu Y Zr Nb Mo Tc Ru Rh Pd Ag La Hf Ta W Re Os Ir Pt Au ',
    'Lanl2-[10s8p7d3f2g]': '  Rh ',
    'Lanl2-[5s4p4d2f]': '  Rh ',
    'Lanl2-[6s4p4d2f]': '  Rh ',
    'Lanl2DZ+1d1f': '  Rh ',
    'Lanl2DZ+2s2p2d2f': '  Rh ',
    'LANL2TZ':
    '  Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd La Hf Ta W Re Os Ir Pt Au Hg ',
    'LANL2TZ+': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'LANL2TZ(f)':
    '  Sc Ti V Cr Mn Fe Co Ni Cu Y Zr Nb Mo Tc Ru Rh Pd Ag La Hf Ta W Re Os Ir Pt Au ',
    'm6-31G': '  Sc Ti V Cr Mn Fe Co Ni Cu ',
    'm6-31G*': '  Sc Ti V Cr Mn Fe Co Ni Cu ',
    'maug-cc-pV(D+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'maug-cc-pVDZ': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'maug-cc-pV(Q+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'maug-cc-pVQZ': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'maug-cc-pV(T+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'maug-cc-pVTZ': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'may-cc-pV(Q+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'may-cc-pV(T+d)Z': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'McLean/Chandler VTZ': '  Na Mg Al Si P S Cl Ar ',
    'MG3S': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'MIDI!': '  H C N O F Si P S Cl Br I ',
    'MIDI(Huzinaga)': '  H He Li Be B C N O F Ne Na Al Si P S Cl Ar K Cs ',
    'MINI(Huzinaga)': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca ',
    'MINI(Scaled)': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca ',
    'modified LANL2DZ':
    '  Sc Ti V Cr Mn Fe Co Ni Cu Y Zr Nb Mo Tc Ru Rh Pd Ag La Hf Ta W Re Os Ir Pt Au ',
    'NASA Ames ANO': '  H B C N O F Ne Al P Ti Fe Ni ',
    'NASA Ames ANO2': '  Sc Ti V ',
    'NASA Ames cc-pCV5Z': '  Ti ',
    'NASA Ames cc-pCVQZ': '  Ti ',
    'NASA Ames cc-pCVTZ': '  Ti ',
    'NASA Ames cc-pV5Z': '  Ti Fe ',
    'NASA Ames cc-pVQZ': '  Ti Fe ',
    'NASA Ames cc-pVTZ': '  Ti Fe ',
    'Partridge Uncontracted 1':
    '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr ',
    'Partridge Uncontracted 2':
    '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'Partridge Uncontracted 3':
    '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'Partridge Uncontracted 4': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'pc-0':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Ga Ge As Se Br Kr ',
    'pc-1':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'pc-2':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'pc-3':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'pc-4':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'pcemd-2': '  He H Be Li Ne C N F B O Mg Na Ar P Cl Al Si S ',
    'pcemd-3': '  He H Be Li Ne C N F B O Mg Na Ar P Cl Al Si S ',
    'pcemd-4': '  He H Be Li Ne C N F B O Mg Na Ar P Cl Al Si S ',
    'pcJ-0': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'pcJ-0_2006': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'pcJ-1': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'pcJ-1_2006': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'pcJ-2': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'pcJ-2_2006': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'pcJ-3': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'pcJ-3_2006': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'pcJ-4': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'pcJ-4_2006': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'pcS-0': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'pcS-1': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'pcS-2': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'pcS-3': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'pcS-4': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'pSBKJC': '  C N O F Si P S Cl Ge As Se Br Sn Sb Te I ',
    'Pt': '  Pt ',
    'pV6Z': '  H C N O ',
    'pV7Z': '  H C N O F Ne S ',
    'Roos_ANO_DZ':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'Roos_ANO_TZ':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'Roos Augmented Double Zeta ANO':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'Roos Augmented Triple Zeta ANO':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    's3-21G': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    's3-21G*': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    's6-31G': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    's6-31G*': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'Sadlej pVTZ': '  H Li Be C N O F Na Mg Si P S Cl K Ca Br Rb Sr I ',
    'SDB-aug-cc-pVQZ': '  Ga Ge As Se Br In Sn Sb Te I ',
    'SDB-aug-cc-pVTZ': '  Ga Ge As Se Br In Sn Sb Te I ',
    'SDB-cc-pVQZ': '  Ga Ge As Se Br Kr In Sn Sb Te I Xe ',
    'SDB-cc-pVTZ': '  Ga Ge As Se Br Kr In Sn Sb Te I Xe ',
    'STO-2G': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sr ',
    'STO-3G':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I ',
    'STO-3G*': '  Na Mg Al Si P S Cl Ar ',
    'STO-6G':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'SV+Double': '  B C N O F Ne ',
    'SV+Rydberg': '  B C N O F Ne ',
    'SV': '  H Li Be B C N O F Ne ',
    'SVP+Rydberg': '  B C N O F Ne ',
    'SVP': '  H Li Be B C N O F Ne ',
    'TZ': '  H Li Be B C N O F Ne ',
    'TZVP': '  H C N O F Al Si P S Cl Ar ',
    'UGBS':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pu Am Cf Es Fm Md No Lr ',
    'un-ccemd-ref': '  He H Be Li Ne C N F B O Mg Na Ar P Cl Al Si S ',
    'un-pcemd-ref': '  He H Be Li Ne C N F B O Mg Na Ar P Cl Al Si S ',
    'Wachters+f': '  Sc Ti V Cr Mn Fe Co Ni Cu ',
    'WTBS':
    '  He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Z3Pol': '  H C N O F Si P S Cl ',
    'aug-cc-pV5Z-PP':
    '  Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'aug-cc-pVDZ-PP':
    '  Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'aug-cc-pVQZ-PP':
    '  Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'aug-cc-pVTZ-PP':
    '  Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'cc-pV5Z-PP':
    '  Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'cc-pVDZ-PP':
    '  Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'cc-pVQZ-PP':
    '  Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'cc-pVTZ-PP':
    '  Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'cc-pwCV5Z-PP':
    '  Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg Ga Ge As Se Br Kr In Sn Sb Te I Xe Tl Pb Bi Po At Rn ',
    'cc-pwCVDZ-PP':
    '  Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'cc-pwCVQZ-PP':
    '  Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'cc-pwCVTZ-PP':
    '  Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'CRENBL_ECP':
    '  H Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Uun Uuu Uub Uut Uuq Uup Uuh Uus ',
    'CRENBS ECP':
    '  Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Rf Db Sg Bh Hs Mt Uun Uuu Uub Uut Uuq Uup Uuh Uus ',
    'Def2-QZVPP':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Def2-QZVPPD':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Def2-TZVPP':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Def2-TZVPPD':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'dhf-QZVPP':
    '  Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'dhf-TZVPP':
    '  Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Hay-Wadt MB (n+1) ECP':
    '  K Ca Sc Ti V Cr Mn Fe Co Ni Cu Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cs Ba La Hf Ta W Re Os Ir Pt Au ',
    'Hay-Wadt VDZ (n+1) ECP':
    '  K Ca Sc Ti V Cr Mn Fe Co Ni Cu Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cs Ba La Ta W Re Os Ir Pt Au ',
    'LANL2DZ ECP':
    '  H Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Pb Bi U Np Pu ',
    'LANL2DZdp ECP': '  H C N O F Si P S Cl Ge As Se Br Sn Sb Te I Pb Bi ',
    'SBKJC VDZ ECP':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'Stuttgart RLC ECP':
    '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Zn Ga Ge As Se Br Kr Rb Sr In Sn Sb Te I Xe Cs Ba Hg Tl Pb Bi Po At Rn Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr ',
    'Stuttgart RSC 1997 ECP':
    '  K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Cs Ba Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Hf Ta W Re Os Ir Pt Au Hg Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Db ',
    'Stuttgart RSC ANO/ECP': '  La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu ',
    'Stuttgart RSC Segmented/ECP':
    '  La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu ',
    'aug-cc-pV5Z-DK':
    '  H He B C N O F Ne Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-cc-pVDZ-DK':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-cc-pVQZ-DK':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-cc-pVTZ-DK':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd ',
    'aug-cc-pwCV5Z-DK': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'aug-cc-pwCVQZ-DK': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'aug-cc-pwCVTZ-DK':
    '  Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd ',
    'cc-pV5Z-DK':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'cc-pVDZ-DK':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'cc-pVQZ-DK':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'cc-pV(T+d)Z-DK': '  Al Si P S Cl Ar ',
    'cc-pVTZ-DK':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd ',
    'cc-pwCV5Z-DK': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'cc-pwCVQZ-DK': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'cc-pwCVTZ-DK': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd ',
    'Cologne_DKH2': '  La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu ',
    'SARC-DKH':
    '  La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr ',
    'SARC-ZORA':
    '  La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr ',
    'cc-pV5Z(fi/sf/fw)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pV5Z(fi/sf/lc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pV5Z(fi/sf/sc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pV5Z(pt/sf/fw)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pV5Z(pt/sf/lc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pV5Z(pt/sf/sc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVDZ(fi/sf/fw)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVDZ(fi/sf/lc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVDZ(fi/sf/sc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVDZ(pt/sf/fw)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVDZ(pt/sf/lc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVDZ(pt/sf/sc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVQZ(fi/sf/fw)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVQZ(fi/sf/lc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVQZ(fi/sf/sc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVQZ(pt/sf/fw)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVQZ(pt/sf/lc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVQZ(pt/sf/sc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVTZ(fi/sf/fw)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVTZ(fi/sf/lc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVTZ(fi/sf/sc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVTZ(pt/sf/fw)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVTZ(pt/sf/lc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pVTZ(pt/sf/sc)':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'aug-cc-pV5Z-PP_OPTRI': '  Cu Zn Ag Cd Au Hg ',
    'aug-cc-pV5Z-RI diffuse': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-cc-pV5Z_OPTRI': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-cc-pV6Z-RI diffuse': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-cc-pVDZ-PP_OPTRI': '  Cu Zn Ag Cd Au Hg ',
    'aug-cc-pVDZ-RI diffuse':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'aug-cc-pVDZ_OPTRI': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-cc-pVQZ-PP_OPTRI': '  Cu Zn Ag Cd Au Hg ',
    'aug-cc-pVQZ-RI diffuse':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'aug-cc-pVQZ_OPTRI': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-cc-pVTZ-PP_OPTRI': '  Cu Zn Ag Cd Au Hg ',
    'aug-cc-pVTZ-RI diffuse':
    '  H He B C N O F Ne Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'aug-cc-pVTZ_OPTRI': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'aug-cc-pwCV5Z-PP_OPTRI': '  Cu Zn Ag Cd Au Hg ',
    'aug-cc-pwCVDZ-PP_OPTRI': '  Cu Zn Ag Cd Au Hg ',
    'aug-cc-pwCVQZ-PP_OPTRI': '  Cu Zn Ag Cd Au Hg ',
    'aug-cc-pwCVTZ-PP_OPTRI': '  Cu Zn Ag Cd Au Hg ',
    'cc-pV(5+d)Z-RI': '  Al Si P S Cl Ar ',
    'cc-pV5Z-PP-RI':
    '  Cu Zn Ga Ge As Se Br Kr Ag Cd In Sn Sb Te I Xe Au Hg Tl Pb Bi Po At Rn ',
    'cc-pV5Z-RI': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar ',
    'cc-pV(6+d)Z-RI': '  Al Si P S Cl Ar ',
    'cc-pV6Z-RI': '  H He B C N O F Ne Al Si P S Cl Ar ',
    'cc-pV(D+d)Z-RI': '  Al Si P S Cl Ar ',
    'cc-pVDZ-PP-RI':
    '  Cu Zn Ga Ge As Se Br Kr Ag Cd In Sn Sb Te I Xe Au Hg Tl Pb Bi Po At Rn ',
    'cc-pVDZ-RI':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pV(Q+d)Z-RI': '  Al Si P S Cl Ar ',
    'cc-pVQZ-PP-RI':
    '  Cu Zn Ga Ge As Se Br Kr Ag Cd In Sn Sb Te I Xe Au Hg Tl Pb Bi Po At Rn ',
    'cc-pVQZ-RI':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pV(T+d)Z-RI': '  Al Si P S Cl Ar ',
    'cc-pVTZ-PP-RI':
    '  Cu Zn Ga Ge As Se Br Kr Ag Cd In Sn Sb Te I Xe Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'cc-pVTZ-RI':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ga Ge As Se Br Kr ',
    'cc-pwCV5Z-RI': '  Cu Zn Ag Cd Au Hg ',
    'cc-pwCV5Z-RI tight': '  B C N O F Ne Al Si P S Cl Ar ',
    'cc-pwCVDZ-RI': '  Cu Zn Ag Cd Au Hg ',
    'cc-pwCVDZ-RI tight': '  B C N O F Ne Al Si P S Cl Ar ',
    'cc-pwCVQZ-RI': '  Cu Zn Ag Cd Au Hg ',
    'cc-pwCVQZ-RI tight': '  B C N O F Ne Al Si P S Cl Ar ',
    'cc-pwCVTZ-RI': '  Cu Zn Ag Cd Au Hg ',
    'cc-pwCVTZ-RI tight': '  B C N O F Ne Al Si P S Cl Ar ',
    'Ahlrichs Coulomb Fitting':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At ',
    'aug-cc-pV5Z-PP MP2 Fitting':
    '  Y Zr Nb Mo Tc Ru Rh Pd Hf Ta W Re Os Ir Pt ',
    'aug-cc-pVDZ-PP MP2 Fitting':
    '  Y Zr Nb Mo Tc Ru Rh Pd Hf Ta W Re Os Ir Pt ',
    'aug-cc-pVQZ-PP MP2 Fitting':
    '  Y Zr Nb Mo Tc Ru Rh Pd Hf Ta W Re Os Ir Pt ',
    'aug-cc-pVTZ-PP MP2 Fitting':
    '  Y Zr Nb Mo Tc Ru Rh Pd Hf Ta W Re Os Ir Pt ',
    'aug-cc-pVTZ MP2 Fitting': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'cc-pV5Z-PP MP2 Fitting': '  Y Zr Nb Mo Tc Ru Rh Pd Hf Ta W Re Os Ir Pt ',
    'cc-pVDZ-PP MP2 Fitting': '  Y Zr Nb Mo Tc Ru Rh Pd Hf Ta W Re Os Ir Pt ',
    'cc-pVQZ-PP MP2 Fitting': '  Y Zr Nb Mo Tc Ru Rh Pd Hf Ta W Re Os Ir Pt ',
    'cc-pVTZ-PP MP2 Fitting': '  Y Zr Nb Mo Tc Ru Rh Pd Hf Ta W Re Os Ir Pt ',
    'cc-pVTZ MP2 Fitting': '  Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'cc-pwCV5Z-PP MP2 Fitting': '  Hf Ta W Re Os Ir Pt ',
    'cc-pwCVDZ-PP MP2 Fitting': '  Hf Ta W Re Os Ir Pt ',
    'cc-pwCVQZ-PP MP2 Fitting': '  Hf Ta W Re Os Ir Pt ',
    'cc-pwCVTZ-PP MP2 Fitting': '  Hf Ta W Re Os Ir Pt ',
    'DeMon Coulomb Fitting':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe ',
    'DGauss A1 DFT Coulomb Fitting':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe ',
    'DGauss A1 DFT Exchange_Fitting':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe ',
    'DGauss A2 DFT Coulomb Fitting':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'DGauss A2 DFT Exchange Fitting':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn ',
    'Weigend Coulomb Fitting':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn ',
    'cc-pCV5Z': '  B C N O F Ne Na Mg Al Si P S Cl Ar Ca ',
    'cc-pCV6Z': '  B C N O F Ne Al Si P S Cl Ar ',
    'cc-pCV6Z(old)': '  O ',
    'cc-pCVDZ': '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ca ',
    'cc-pCVDZ(old)': '  Al Si P S Cl Ar ',
    'cc-pCVQZ': '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ca ',
    'cc-pCVQZ(old)': '  Al Si P S Cl Ar ',
    'cc-pCVTZ': '  Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Ca ',
    'cc-pCVTZ(old)': '  Al Si P S Cl Ar ',
    'aug-cc-pV5Z-DK Diffuse':
    '  H He B C N O F Ne Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-cc-pVDZ-DK Diffuse':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-cc-pVQZ-DK Diffuse':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ',
    'aug-cc-pVTZ-DK Diffuse':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Y Zr Nb Mo Tc Ru Rh Pd ',
    'cc-pwCV5Z Tight': '  B C N O F Ne Al Si P S Cl Ar ',
    'cc-pwCVDZ Tight': '  B C N O F Ne Al Si P S Cl Ar ',
    'cc-pwCVQZ Tight': '  B C N O F Ne Al Si P S Cl Ar Br ',
    'cc-pwCVTZ Tight': '  B C N O F Ne Al Si P S Cl Ar '
}

# This list contains all correlation-consistent basis sets, which were designed using
# spherical harmonics. To use them, the spherical keyword should be present in the BASIS directive.
CC_basis_sets = [
    'apr-cc-pV(Q+d)Z', 'aug-cc-pCV5Z', 'aug-cc-pCVDZ', 'aug-cc-pCVQZ',
    'aug-cc-pCV(T+d)Z', 'aug-cc-pCVTZ', 'aug-cc-pV(5+d)Z', 'aug-cc-pV5Z',
    'aug-cc-pV(6+d)Z', 'aug-cc-pV6Z', 'aug-cc-pV(D+d)Z', 'aug-cc-pVDZ',
    'aug-cc-pV(Q+d)Z', 'aug-cc-pVQZ', 'aug-cc-pV(T+d)Z', 'aug-cc-pVTZ',
    'aug-cc-pVTZ-J', 'aug-cc-pwCV5Z', 'aug-cc-pwCV5Z-NR', 'aug-cc-pwCVDZ',
    'aug-cc-pwCVQZ', 'aug-cc-pwCVQZ-NR', 'aug-cc-pwCVTZ', 'aug-cc-pwCVTZ-NR',
    'cc-pCV5Z', 'cc-pCV6Z', 'cc-pCV6Z(old)', 'cc-pCVDZ', 'cc-pCVDZ(old)',
    'cc-pCVQZ', 'cc-pCVQZ(old)', 'cc-pCVTZ', 'cc-pCVTZ(old)', 'cc-pV(5+d)Z',
    'cc-pV5Z', 'cc-pV(6+d)Z', 'cc-pV6Z', 'cc-pV8Z', 'cc-pV9Z', 'cc-pV(D+d)Z',
    'cc-pVDZ', 'cc-pVDZ(seg-opt)', 'cc-pV(Q+d)Z', 'cc-pVQZ',
    'cc-pVQZ(seg-opt)', 'cc-pV(T+d)Z', 'cc-pVTZ', 'cc-pVTZ(seg-opt)',
    'cc-pwCV5Z', 'cc-pwCV5Z', 'cc-pwCV5Z-NR', 'cc-pwCVDZ', 'cc-pwCVQZ',
    'cc-pwCVQZ-NR', 'cc-pwCVTZ', 'cc-pwCVTZ-NR', 'd-aug-cc-pV5Z',
    'd-aug-cc-pV6Z', 'd-aug-cc-pVDZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pVTZ',
    'jul-cc-pV(D+d)Z', 'jul-cc-pV(Q+d)Z', 'jul-cc-pV(T+d)Z', 'jun-cc-pV(D+d)Z',
    'jun-cc-pV(Q+d)Z', 'jun-cc-pV(T+d)Z', 'maug-cc-pV(D+d)Z', 'maug-cc-pVDZ',
    'maug-cc-pV(Q+d)Z', 'maug-cc-pVQZ', 'maug-cc-pV(T+d)Z', 'maug-cc-pVTZ',
    'may-cc-pV(Q+d)Z', 'may-cc-pV(T+d)Z', 'aug-cc-pV5Z-PP', 'aug-cc-pVDZ-PP',
    'aug-cc-pVQZ-PP', 'aug-cc-pVTZ-PP', 'cc-pV5Z-PP', 'cc-pVDZ-PP',
    'cc-pVQZ-PP', 'cc-pVTZ-PP', 'cc-pwCV5Z-PP', 'cc-pwCVDZ-PP', 'cc-pwCVQZ-PP',
    'cc-pwCVTZ-PP', 'aug-cc-pV5Z-DK', 'aug-cc-pVDZ-DK', 'aug-cc-pVQZ-DK',
    'aug-cc-pVTZ-DK', 'aug-cc-pwCV5Z-DK', 'aug-cc-pwCVQZ-DK',
    'aug-cc-pwCVTZ-DK', 'cc-pV5Z-DK', 'cc-pVDZ-DK', 'cc-pVQZ-DK',
    'cc-pV(T+d)Z-DK', 'cc-pVTZ-DK', 'cc-pwCV5Z-DK', 'cc-pwCVQZ-DK',
    'cc-pwCVTZ-DK', 'cc-pV5Z(fi/sf/fw)', 'cc-pV5Z(fi/sf/lc)',
    'cc-pV5Z(fi/sf/sc)', 'cc-pV5Z(pt/sf/fw)', 'cc-pV5Z(pt/sf/lc)',
    'cc-pV5Z(pt/sf/sc)', 'cc-pVDZ(fi/sf/fw)', 'cc-pVDZ(fi/sf/lc)',
    'cc-pVDZ(fi/sf/sc)', 'cc-pVDZ(pt/sf/fw)', 'cc-pVDZ(pt/sf/lc)',
    'cc-pVDZ(pt/sf/sc)', 'cc-pVQZ(fi/sf/fw)', 'cc-pVQZ(fi/sf/lc)',
    'cc-pVQZ(fi/sf/sc)', 'cc-pVQZ(pt/sf/fw)', 'cc-pVQZ(pt/sf/lc)',
    'cc-pVQZ(pt/sf/sc)', 'cc-pVTZ(fi/sf/fw)', 'cc-pVTZ(fi/sf/lc)',
    'cc-pVTZ(fi/sf/sc)', 'cc-pVTZ(pt/sf/fw)', 'cc-pVTZ(pt/sf/lc)',
    'cc-pVTZ(pt/sf/sc)', 'aug-cc-pV5Z-PP_OPTRI', 'aug-cc-pV5Z-RI diffuse',
    'aug-cc-pV5Z_OPTRI', 'aug-cc-pV6Z-RI diffuse', 'aug-cc-pVDZ-PP_OPTRI',
    'aug-cc-pVDZ-RI diffuse', 'aug-cc-pVDZ_OPTRI', 'aug-cc-pVQZ-PP_OPTRI',
    'aug-cc-pVQZ-RI diffuse', 'aug-cc-pVQZ_OPTRI', 'aug-cc-pVTZ-PP_OPTRI',
    'aug-cc-pVTZ-RI diffuse', 'aug-cc-pVTZ_OPTRI', 'aug-cc-pwCV5Z-PP_OPTRI',
    'aug-cc-pwCVDZ-PP_OPTRI', 'aug-cc-pwCVQZ-PP_OPTRI',
    'aug-cc-pwCVTZ-PP_OPTRI', 'cc-pV(5+d)Z-RI', 'cc-pV5Z-PP-RI', 'cc-pV5Z-RI',
    'cc-pV(6+d)Z-RI', 'cc-pV6Z-RI', 'cc-pV(D+d)Z-RI', 'cc-pVDZ-PP-RI',
    'cc-pVDZ-RI', 'cc-pV(Q+d)Z-RI', 'cc-pVQZ-PP-RI', 'cc-pVQZ-RI',
    'cc-pV(T+d)Z-RI', 'cc-pVTZ-PP-RI', 'cc-pVTZ-RI', 'cc-pwCV5Z-RI',
    'cc-pwCV5Z-RI tight', 'cc-pwCVDZ-RI', 'cc-pwCVDZ-RI tight', 'cc-pwCVQZ-RI',
    'cc-pwCVQZ-RI tight', 'cc-pwCVTZ-RI', 'cc-pwCVTZ-RI tight',
    'aug-cc-pV5Z-PP MP2 Fitting', 'aug-cc-pVDZ-PP MP2 Fitting',
    'aug-cc-pVQZ-PP MP2 Fitting', 'aug-cc-pVTZ-PP MP2 Fitting',
    'aug-cc-pVTZ MP2 Fitting', 'cc-pV5Z-PP MP2 Fitting',
    'cc-pVDZ-PP MP2 Fitting', 'cc-pVQZ-PP MP2 Fitting',
    'cc-pVTZ-PP MP2 Fitting', 'cc-pVTZ MP2 Fitting',
    'cc-pwCV5Z-PP MP2 Fitting', 'cc-pwCVDZ-PP MP2 Fitting',
    'cc-pwCVQZ-PP MP2 Fitting', 'cc-pwCVTZ-PP MP2 Fitting', 'cc-pCV5Z',
    'cc-pCV6Z', 'cc-pCV6Z(old)', 'cc-pCVDZ', 'cc-pCVDZ(old)', 'cc-pCVQZ',
    'cc-pCVQZ(old)', 'cc-pCVTZ', 'cc-pCVTZ(old)', 'aug-cc-pV5Z-DK Diffuse',
    'aug-cc-pVDZ-DK Diffuse', 'aug-cc-pVQZ-DK Diffuse',
    'aug-cc-pVTZ-DK Diffuse', 'cc-pwCV5Z Tight', 'cc-pwCVDZ Tight',
    'cc-pwCVQZ Tight', 'cc-pwCVTZ Tight'
]

# This list contains the names of all basis set available in NWChemEx
# (in the convention of NWChemEx)
nwchemex_basis_list = [
    '3-21g', '4-31g', '5-21g', '6-21g', '6-31++g', '6-31++g_st_',
    '6-31++g_st__st_', '6-31+g', '6-31+g_st_', '6-31+g_st__st_',
    '6-311++g(2d,2p)', '6-311++g(3df,3pd)', '6-311++g', '6-311++g_st_',
    '6-311++g_st__st_', '6-311+g(2d,p)', '6-311+g', '6-311+g_st_',
    '6-311+g_st__st_', '6-311g(2df,2pd)', '6-311g(d,p)', '6-311g',
    '6-311g_st_', '6-311g_st__st_', '6-31g(2df,p)', '6-31g(3df,3pd)',
    '6-31g(d,p)', '6-31g', '6-31g_st_', '6-31g_st__st_', 'aug-cc-pcv5z',
    'aug-cc-pcvdz', 'aug-cc-pcvqz', 'aug-cc-pcvtz', 'aug-cc-pv(5+d)z',
    'aug-cc-pv(d+d)z', 'aug-cc-pv(q+d)z', 'aug-cc-pv(t+d)z',
    'aug-cc-pv5z-optri', 'aug-cc-pv5z-pp-optri', 'aug-cc-pv5z-pp-rifit',
    'aug-cc-pv5z-rifit', 'aug-cc-pv5z', 'aug-cc-pv6z-rifit', 'aug-cc-pv6z',
    'aug-cc-pv7z', 'aug-cc-pvdz-optri', 'aug-cc-pvdz-pp-optri',
    'aug-cc-pvdz-pp-rifit', 'aug-cc-pvdz-rifit', 'aug-cc-pvdz',
    'aug-cc-pvqz-optri', 'aug-cc-pvqz-pp-optri', 'aug-cc-pvqz-pp-rifit',
    'aug-cc-pvqz-rifit', 'aug-cc-pvqz', 'aug-cc-pvtz-optri',
    'aug-cc-pvtz-pp-optri', 'aug-cc-pvtz-pp-rifit', 'aug-cc-pvtz-rifit',
    'aug-cc-pvtz', 'aug-cc-pwcv5z-pp-optri', 'aug-cc-pwcv5z-pp-rifit',
    'aug-cc-pwcv5z-rifit', 'aug-cc-pwcv5z', 'aug-cc-pwcvdz-pp-optri',
    'aug-cc-pwcvdz-pp-rifit', 'aug-cc-pwcvdz-rifit', 'aug-cc-pwcvdz',
    'aug-cc-pwcvqz-pp-optri', 'aug-cc-pwcvqz-pp-rifit', 'aug-cc-pwcvqz-rifit',
    'aug-cc-pwcvqz', 'aug-cc-pwcvtz-pp-optri', 'aug-cc-pwcvtz-pp-rifit',
    'aug-cc-pwcvtz-rifit', 'aug-cc-pwcvtz', 'aug-pv7z', 'cc-pcv5z', 'cc-pcvdz',
    'cc-pcvqz', 'cc-pcvtz', 'cc-pv(5+d)z', 'cc-pv(d+d)z', 'cc-pv(q+d)z',
    'cc-pv(t+d)z', 'cc-pv5z-jkfit', 'cc-pv5z-pp-rifit', 'cc-pv5z-rifit',
    'cc-pv5z', 'cc-pv6z-rifit', 'cc-pv6z', 'cc-pv8z', 'cc-pv9z',
    'cc-pvdz(seg-opt)', 'cc-pvdz-pp-rifit', 'cc-pvdz-rifit', 'cc-pvdz',
    'cc-pvqz(seg-opt)', 'cc-pvqz-jkfit', 'cc-pvqz-pp-rifit', 'cc-pvqz-rifit',
    'cc-pvqz', 'cc-pvtz(seg-opt)', 'cc-pvtz-jkfit', 'cc-pvtz-pp-rifit',
    'cc-pvtz-rifit', 'cc-pvtz', 'cc-pwcv5z-pp-rifit', 'cc-pwcv5z-rifit',
    'cc-pwcv5z', 'cc-pwcvdz-pp-rifit', 'cc-pwcvdz-rifit', 'cc-pwcvdz',
    'cc-pwcvqz-pp-rifit', 'cc-pwcvqz-rifit', 'cc-pwcvqz', 'cc-pwcvtz-pp-rifit',
    'cc-pwcvtz-rifit', 'cc-pwcvtz', 'd-aug-cc-pv5z', 'd-aug-cc-pv6z',
    'd-aug-cc-pvdz', 'd-aug-cc-pvqz', 'd-aug-cc-pvtz', 'pv6z', 'pv7z',
    'sto-2g', 'sto-3g', 'sto-3g_st_', 'sto-4g', 'sto-5g', 'sto-6g'
]

small_basis_elements = {
    '3-21G':
    '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs ',
    '6-31++G*': '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca '
}

#for basis_name in basis_elements:
for basis_name in small_basis_elements:

    nwchemex_basis_name = basis_name.lower().replace("*",
                                                     "_st_").replace(" ", "_")

    if nwchemex_basis_name in nwchemex_basis_list:

        print("Generating atomic density matrices for basis set",
              nwchemex_basis_name, ":")
        outfile = nwchemex_basis_name.replace("_st_", "_star") + ".dat"

        # remove the old files if exist
        if os.path.isfile(outfile):
            os.remove(outfile)

        str1 = (basis_elements[basis_name].strip()).split()
        for cc in str1:
            with open('tmp.nw', 'w') as nwInp:
                if (atom_sym2num(cc) % 2 == 0):
                    nwInp.write("echo\n")
                    nwInp.write("start test-scf-sad\n")
                    nwInp.write(
                        "title \"Test atom HF calculation in {} to check SAD guess\"\n"
                        .format(basis_name))
                    nwInp.write("\n")
                    nwInp.write("scratch_dir ./scratch\n")
                    nwInp.write("permanent_dir ./perm\n")
                    nwInp.write("\n")
                    nwInp.write("charge 0\n")
                    nwInp.write("\n")
                    nwInp.write("geometry\n")
                    nwInp.write(
                        "  {}      0.00000000    0.00000000    0.00000000\n".
                        format(cc))
                    nwInp.write("end\n")
                    nwInp.write("\n")
                    if basis_name in CC_basis_sets:
                        nwInp.write("basis spherical\n")
                    else:
                        nwInp.write("basis\n")
                    nwInp.write("  * library {} \n".format(basis_name))
                    nwInp.write("end\n")
                    nwInp.write("\n")
                    nwInp.write("scf\n")
                    nwInp.write(" singlet\n")
                    nwInp.write(" noscf\n")
                    nwInp.write(" print \"atomic guess density\"\n")
                    nwInp.write("end\n")
                    nwInp.write("\n")
                    nwInp.write("task scf energy\n")
                else:
                    nwInp.write("echo\n")
                    nwInp.write("start test-scf-sad\n")
                    nwInp.write(
                        "title \"Test atom HF calculation in {} to check SAD guess\"\n"
                        .format(basis_name))
                    nwInp.write("\n")
                    nwInp.write("scratch_dir ./scratch\n")
                    nwInp.write("permanent_dir ./perm\n")
                    nwInp.write("\n")
                    nwInp.write("charge 0\n")
                    nwInp.write("\n")
                    nwInp.write("geometry\n")
                    nwInp.write(
                        "  {}      0.00000000    0.00000000    0.00000000\n".
                        format(cc))
                    nwInp.write("end\n")
                    nwInp.write("\n")
                    if basis_name in CC_basis_sets:
                        nwInp.write("basis spherical\n")
                    else:
                        nwInp.write("basis\n")
                    nwInp.write("  * library {} \n".format(basis_name))
                    nwInp.write("end\n")
                    nwInp.write("\n")
                    nwInp.write("scf\n")
                    nwInp.write(" doublet\n")
                    nwInp.write(" noscf\n")
                    nwInp.write(" print \"atomic guess density\"\n")
                    nwInp.write("end\n")
                    nwInp.write("\n")
                    nwInp.write("task scf energy\n")
            subprocess.run(["./mynw.sh", 'tmp'])
            with open('tmp.nwo', "r") as nwOut:
                for line in nwOut:
                    if line.startswith(" global array: Guess"):
                        dim = int((line.replace(',', ":")).split(":")[2])
                        nwOut.readline()
                        nwOut.readline()
                        nwOut.readline()
                        dm_vec = np.zeros(dim * dim)
                        n1 = dim // 6
                        r1 = dim % 6
                        if n1 == 0:
                            for i in range(dim):
                                line1 = ((nwOut.readline()).split())[1:]
                                for j in range(dim):
                                    dm_vec[i * dim + j] = float(line1[j])
                                i += 1
                        else:
                            nn = 0
                            while nn < n1:
                                for i in range(dim):
                                    line1 = ((nwOut.readline()).split())[1:]
                                    for j in range(6):
                                        dm_vec[i * dim + j + 6 * nn] = float(
                                            line1[j])
                                nn += 1
                                nwOut.readline()
                                nwOut.readline()
                                nwOut.readline()
                            # The last partial block of the matrtix
                            if r1 != 0:
                                for i in range(dim):
                                    line1 = (nwOut.readline()).split()[1:]
                                    for j in range(r1):
                                        dm_vec[i * dim + j + 6 * n1] = float(
                                            line1[j])
                                    i += 1
                        with open(outfile, "a") as dmOut:
                            dmOut.write("{}\n".format(cc))
                            ii = 0
                            for el in dm_vec:
                                dmOut.write("{:12.8f} ".format(el))
                                ii += 1
                                if ii % 6 == 0:
                                    dmOut.write("\n")
                            dmOut.write("\n\n")
