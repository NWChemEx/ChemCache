#!/usr/bin/python3

""" Python script to generate all atomic density matrices for a specifice general basis set 
    (not in the NWChem library) with NWChem calculations

  To properly run this script, one need to have a NWChem basis set file 
  (setting below: nw_basis_file = basis_name+".nw") under the working directory. Usually this
  file can be obtained from the basis set exchange website.
  In addition, one needs to set a dictionary containing the basis set name and all elements
  associated to it (setting below: small_basis_elements = { *** }.
  Once all those settings are done, simply run the script as "./genbase_nwchem_dm.py" in the
  command line.
#
# The script write a NWChem input file (tmp.nw) for each atom in the value list corresponding to
# the basis name key, and call NWChem to run a simple SAD guess generation calculaton (noscf).
#
# After the calculation the SAD guess (atomic density matrix) is extracted and reformmated, and
# stored accumulatively in a "basis_name.dat" file (e. g., STO-3G.dat).
#
# The atomic density matrices are stored as text blocks representing the 1-D flatten version of
# the corresponding 2-D matrices.
"""


import numpy as np
import subprocess
import os

def atom_sym2num(sym):
    elements = {'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,
            'Ne':10,'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,
            'Ar':18,'K':19,'Ca':20,'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,
            'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,'Ga':31,'Ge':32,'As':33,
            'Se':34,'Br':35,'Kr':36,'Rb':37,'Sr':38,'Y':39,'Zr':40,'Nb':41,
            'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,
            'Sn':50,'Sb':51,'Te':52,'I':53,'Xe':54,'Cs':55,'Ba':56,'La':57,
            'Ce':58,'Pr':59,'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,
            'Dy':66,'Ho':67,'Er':68,'Tm':69,'Yb':70,'Lu':71,'Hf':72,'Ta':73,
            'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,'Tl':81,
            'Pb':82,'Bi':83,'Po':84,'At':85,'Rn':86,'Fr':87,'Ra':88,'Ac':89,
            'Th':90,'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,
            'Cf':98,'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,'Rf':104,
            'Db':105,'Sg':106,'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,
            'Cn':112,'Nh':113,'Fl':114,'Mc':115,'Lv':116,'Ts':117,'Og':118}

    try:
        return elements[sym]
    except:
        print("Cannot find the atomic number for",sym,"!")

# To get the atomic symbol from the corresponding atomic number
#print(list(elements.keys())[list(elements.values()).index(3)])

#small_basis_elements = {'mini' : '  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn'}

#small_basis_elements = {'mini' : '  H He Li'}
small_basis_elements = {'mini' : '  La'}

for basis_name in small_basis_elements:

    #nwchemex_basis_name = basis_name.lower().replace("*","_st_").replace(" ","_")


    print("Generating atomic density matrices for basis set",basis_name,":")
    outfile = basis_name+".dat"

    # basis set file (in nwchem format)
    nw_basis_file = basis_name+".nw" # hardwired here
    
    # remove the old files if exist
    if os.path.isfile(outfile):
        os.remove(outfile)
    
    str1 = (small_basis_elements[basis_name].strip()).split()
    for cc in str1:
        print ("Running NWChem calculation for atom "+cc)
        with open('tmp.nw', 'w') as nwInp:
            if (atom_sym2num(cc)%2 == 0):
                nwInp.write("echo\n")
                nwInp.write("start test-scf-sad\n")
                nwInp.write("title \"Test atom HF calculation in {} to check SAD guess\"\n".format(basis_name))
                nwInp.write("\n")
                nwInp.write("scratch_dir ./scratch\n")
                nwInp.write("permanent_dir ./perm\n")
                nwInp.write("\n")
                nwInp.write("charge 0\n")
                nwInp.write("\n")
                nwInp.write("geometry\n")
                nwInp.write("  {}      0.00000000    0.00000000    0.00000000\n".format(cc))
                nwInp.write("end\n")
                nwInp.write("\n")
                # Read the general basis set from external basis set file
                with open(nw_basis_file,"r") as nwBasis:
                    inBlock = 0
                    for line in nwBasis:
                        line1 = line.split()
                        if line[0] == "#" and line[1] != "B":
                                nwInp.write(line)
                        elif line == "\n" or line[0:5] == "BASIS":
                            nwInp.write(line)
                        elif (line1[0] == cc or inBlock == 1) and (line[0:6] != "#BASIS"):
                            inBlock = 1
                            nwInp.write(line)
                        elif inBlock == 1 and line[0:6] == "#BASIS":
                            nwInp.write("END\n")
                            break
                        else:
                            pass
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
                nwInp.write("title \"Test atom HF calculation in {} to check SAD guess\"\n".format(basis_name))
                nwInp.write("\n")
                nwInp.write("scratch_dir ./scratch\n")
                nwInp.write("permanent_dir ./perm\n")
                nwInp.write("\n")
                nwInp.write("charge 0\n")
                nwInp.write("\n")
                nwInp.write("geometry\n")
                nwInp.write("  {}      0.00000000    0.00000000    0.00000000\n".format(cc))
                nwInp.write("end\n")
                nwInp.write("\n")
                # Read the general basis set from external basis set file
                with open(nw_basis_file,"r") as nwBasis:
                    inBlock = 0
                    for line in nwBasis:
                        line1 = line.split()
                        if line[0] == "#" and line[1] != "B":
                                nwInp.write(line)
                        elif line == "\n" or line[0:5] == "BASIS":
                            nwInp.write(line)
                        elif (line1[0] == cc or inBlock == 1) and (line[0:6] != "#BASIS"):
                            inBlock = 1
                            nwInp.write(line)
                        elif inBlock == 1 and line[0:6] == "#BASIS":
                            nwInp.write("END\n")
                            break
                        else:
                            pass
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
                    dim = int((line.replace(',',":")).split(":")[2])
                    nwOut.readline()
                    nwOut.readline()
                    nwOut.readline()
                    dm_vec = np.zeros(dim*dim)
                    n1 = dim//6
                    r1 = dim%6
                    if n1 == 0:
                        for i in range(dim) :
                            line1 = ((nwOut.readline()).split())[1:]
                            for j in range(dim):
                                dm_vec[i*dim+j] = float(line1[j])
                            i += 1
                    else:
                        nn = 0
                        while nn < n1:
                            for i in range(dim):
                                line1 = ((nwOut.readline()).split())[1:]
                                for j in range(6):
                                    dm_vec[i*dim+j+6*nn] = float(line1[j])
                            nn += 1
                            nwOut.readline()
                            nwOut.readline()
                            nwOut.readline()
                        # The last partial block of the matrtix
                        if r1 != 0:
                            for i in range(dim) :
                                line1 = (nwOut.readline()).split()[1:]
                                for j in range(r1):
                                    dm_vec[i*dim+j+6*n1] = float(line1[j])
                                i += 1
                    with open(outfile, "a") as dmOut:
                        dmOut.write("{}\n".format(cc))
                        ii = 0
                        for el in dm_vec:
                            dmOut.write("{:12.8f} ".format(el))
                            ii += 1
                            if ii%6 == 0:
                                dmOut.write("\n")
                        dmOut.write("\n\n")
                                
