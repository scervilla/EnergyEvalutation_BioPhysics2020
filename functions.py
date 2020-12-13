import math
from math import e, sqrt
import matplotlib.pyplot as pl

def vdw_energy(atom1, atom2):
    '''Function that calculates the van der Waals energy between two atoms'''
    sig1 = vars(atom1.xtra['vdw'])['sig']
    sig2 = vars(atom2.xtra['vdw'])['sig']
    dis = atom1-atom2
    return 4*sqrt(vars(atom1.xtra['vdw'])['eps']*vars(atom2.xtra['vdw'])['eps'])*(sig1**6*sig2**6/dis**12 - sig1**3*sig2**3/dis**6)

def elec_energy(atom1, atom2, constant):
    '''Function that calculates the electrostatics energy between two atoms.
       The third parameter specify the dialectric constant (three options: vacuum, water, and Mehler-Solmajer(MS)).'''
    dis = atom1-atom2
    if constant == 'vaccum':
        eps = 1
    elif constant == 'water':
        eps = 80
    else: 
        eps = 86.9525 / (1-7.7839*e**(-0.3153*dis)) - 8.5525
    return 332.16*atom1.xtra['charge']*atom2.xtra['charge']/(dis*eps)

def solv_energy(atom):
    '''Fucntion that calculates the solvation energy through ASA surface of an atom (except for hydrogens)'''
    if atom.xtra['atom_type'] not in ('HD', 'H'):
        return float(atom.xtra['EXP_NACCESS'])*vars(atom.xtra['vdw'])['fsrf']
    else:
        return 0

def create_plot(energies):
    '''Create a representation of ΔΔG for each position of the mutation. It is added a level for each point an horizontal line in y = 0.'''
    change = list()
    x = sorted(energies)
    for pos in x:
        change.append(sum(energies[pos]))    
    pl.plot(x,change)
    for i, txt in enumerate(x):
        pl.annotate(txt, (x[i], change[i]))
    pl.axhline(y=0, color = 'red', linestyle='--')
    pl.xlabel("Mutation position")
    pl.ylabel("ΔΔG")
    pl.show()

def get_unfold_file():
    '''Stores solvation energy of aminoacids in unfolded state from the file unfolded_e.txt'''
    f = open('data/unfolded_e.txt', 'r')
    unfold = dict()
    for line in f:
        aa, val = line.split()
        unfold[aa] = float(val)
    f.close()    
    return unfold


def unfolded_e(st, unfold):
    '''Calculates solvation energy for the original protein'''
    count = 0
    for res in st.get_residues():
        if res.get_resname() == 'HIE':
            count += unfold['HIS']
        elif res.get_resname() != 'ACE' and res.get_resname() != 'NME':
            count += unfold[res.get_resname()]
    return count

def eval_energies(atoms, st, u):
    '''Evaluates Vdw, electrostatics and solvent energies of the orginial protein and each of the position mutation. It returns a tuple with a dictioanry with the position 
    of the mutation as key and a list with elec, VDW, solvent (folded) and solvent(unfolded). The parameters corresponds to: list of atoms, structure of the protein and 
    the dictionary of energies in an unfolded state.''' 
    i = 1
    count_vdw = 0
    count_elec = 0
    count_solv = 0
    energies = dict()
    print('Evaluating energies')
    print('0%')
    
    for res in st.get_residues(): 
        if res.get_resname() != "ALA" and res.get_resname() != "GLY" and res.get_resname() != 'NME' and res.get_resname() != 'ACE':
            if res.get_resname() == 'HIE':
                energies[res.get_id()[1]] = [0, 0, 0, u['ALA']-u['HIS']]
            else:
                energies[res.get_id()[1]] = [0, 0, 0, u['ALA']-u[res.get_resname()]]
    l = len(atoms)
    if l % 2 != 0:
        l = l-1
    for at1 in atoms:
        for at2 in atoms[i:]:
            if at2.get_parent() != at1.get_parent() and not (vars(at1)['name'] == 'C' and vars(at2)['name'] == 'N' and vars(at1)['full_id'][3][1]+1 == vars(at2)['full_id'][3][1]):
                if 2.5  < at1-at2:
                    vdw = vdw_energy(at1, at2)

                    count_vdw += vdw
                    if vars(at1)['name'] not in ('N', 'CA', 'C', 'O', 'CB') and at1.get_parent().get_id()[1] in energies: 
                        energies[at1.get_parent().get_id()[1]][1] -= vdw
                    if vars(at2)['name'] not in ('N', 'CA', 'C', 'O', 'CB') and at2.get_parent().get_id()[1] in energies:
                        energies[at2.get_parent().get_id()[1]][1] -= vdw
                elec = elec_energy(at1, at2, 'MS')
                count_elec += elec 
                if vars(at1)['name'] not in ('N', 'CA', 'C', 'O', 'CB') and at1.get_parent().get_id()[1] in energies:
                    energies[at1.get_parent().get_id()[1]][0] -= elec
                if vars(at1)['name'] not in ('N', 'CA', 'C', 'O', 'CB') and at2.get_parent().get_id()[1] in energies:
                    energies[at2.get_parent().get_id()[1]][0] -= elec
        solv = solv_energy(at1)
        count_solv += solv
        if vars(at1)['name'] not in ('N', 'CA', 'C', 'O', 'CB') and at1.get_parent().get_id()[1] in energies:
            energies[at1.get_parent().get_id()[1]][2] -= solv
        i += 1

        if i / l == 0.5:
            print('50%')
    print('Completed')
    return (energies, count_elec, count_vdw, count_solv)

def energy_file(energies, principal, pdb, st):
    '''Creates a txt file with all energies evaluation of ALA mutation and the original protein. It contains as columns: Residue mutated, position of the residue, Solvation 
    energy(unfolded), Van der Waals energy, Electrostatic energy, Solvation energy (folded) and ΔΔG'''
    resname = list()
    for res in st.get_residues():
        if res.get_resname() != "ALA" and res.get_resname() != "GLY" and res.get_resname() != 'NME' and res.get_resname() != 'ACE' :
            resname.append(res.get_resname())
    f = open('energy_' + pdb[:-3] + 'txt', 'w+')
    f.write('-Energy Analysis for ALA mutations-\n')
    f.write('PDB: ' + pdb[:4] + '\n')
    f.write('Res\tPos\tUnfolded Energy\t    VDW Energy\t        Electrostatics En.\tSolvation Energy\t ΔΔG\n')
    f.write('Ori\tNon\t'+ str(principal[0]) + '\t' + str(principal[1]) + '\t'+ str(principal[2]) + '\t' + str(principal[3]) + '\n')
    i = 0
    for pos in sorted(energies):
        f.write(resname[i] + '\t' + str(pos) + '\t' + str(energies[pos][3]+principal[0]) + '\t' + str(energies[pos][1]+principal[1]) + '\t' + str(energies[pos][0]+principal[2]) 
        + '\t' + str(energies[pos][2]+principal[3])+ '\t' + str(sum(energies[pos])) + '\n')
        i += 1
    f.close()