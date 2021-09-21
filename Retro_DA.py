#-------------------------------------------------------------------------------
# Copyright (c) 2021 CMCDD Research Group, Bienfait Isamura, Kevin Lobb        #
#                    Chemistry Dept, Rhodes University                         #
# Permission is hereby granted, free of charge, to any person obtaining a      #
# copy of this program and associated documentation files (the "Program"),     #
# to deal in the Software without restriction, including without limitation    #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,     #
# and/or sell copies of the Software, and to permit persons to whom the        #
# Software is furnished to do so, subject to the following conditions:         #
#                                                                              #
# The above copyright notice and this permission notice shall be included in   #
# all copies or substantial portions of the Software.                          #
#                                                                              #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS      #
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,  #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE  #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS #
# IN THE SOFTWARE.                                                             #
#-------------------------------------------------------------------------------

from rdkit import Chem
from rdkit.Chem import Draw,RWMol,rdchem,BondType,AllChem
from rdkit.Chem import rdMolDescriptors,rdMolTransforms
import numpy as np
   
    
def identify_cyclohxn_substructures(mol):
    
    cyclo=Chem.MolFromSmarts('C1=CCCCC1')
    
    substructures=mol.GetSubstructMatches(cyclo,useChirality=True)  # Identifying indexes of atoms in the cyclohxn ring
    
    return substructures, len(substructures)
    


def locate_bonds_of_interest(mol):
    
    count_substructures=identify_cyclohxn_substructures(mol)[1]
    
    all_substructures=identify_cyclohxn_substructures(mol)[0]
    
    result=[]
    
    for i in range(count_substructures):

        result.append([])        
        
        ith_substructure=all_substructures[i]
        
        bonds_of_dienophile_to_modify=[]
        
        index_atoms_in_bonds_to_break=[]
        
        index_atoms_in_initial_ene_double_bond=[]
        
        for j in range(len(ith_substructure)-1):
            
            index_atom1=ith_substructure[j]
            
            index_atom2=ith_substructure[j+1]
            
            if mol.GetBondBetweenAtoms(index_atom1,index_atom2).GetBondType()==2:  #I locate the unique double bond in the cyclohexene substructure
                
                bonds_of_dienophile_to_modify.append(mol.GetBondBetweenAtoms(index_atom1,index_atom2))
                
                for atom in [mol.GetAtomWithIdx(x) for x in list([index_atom1,index_atom2])]:
                    
                    for neighbor in atom.GetNeighbors():
                        
                        if (neighbor.GetIdx() in ith_substructure)==True and (neighbor.GetIdx() not in list([index_atom1,index_atom2])):  #Finding the direct next neighbor of one atom in the double bond of the cyclohexene ring. This atom is the begin atom of the bond to break
                            
                            bonds_of_dienophile_to_modify.append(mol.GetBondBetweenAtoms(atom.GetIdx(),neighbor.GetIdx()))
                            
                            index_atoms_in_bonds_to_break.append(neighbor.GetIdx())
                            
                            for next_ring_neighbor in neighbor.GetNeighbors():
                                
                                if (next_ring_neighbor.GetIdx() in ith_substructure)==True and (next_ring_neighbor.GetIdx() not in list([index_atom1,index_atom2])):  #Finding the second next neighbor of one atom in the double bond of the cyclohexene ring at the same side. It is the end atom of the bond to remove
            
                                    index_atoms_in_bonds_to_break.append(next_ring_neighbor.GetIdx())                                                                 # This atom is also an extremity of the single bond which will be turned into double bond after effecting the retro Diels Alder
                                    
                                    index_atoms_in_initial_ene_double_bond.append(next_ring_neighbor.GetIdx())
        
        
        result[i].append(bonds_of_dienophile_to_modify)
        
        result[i].append(index_atoms_in_bonds_to_break)
        
        result[i].append(index_atoms_in_initial_ene_double_bond)
        
    return result
    
    
def process_retro_Diels_Alder(mol):
    
    list_retroDA_adducts=[]
    
    mol1=Chem.MolFromSmiles('C=C-C')
            
    double_bond=mol1.GetBonds()[0]
        
    single_bond=mol1.GetBonds()[1]
    
    for i in range(len(locate_bonds_of_interest(mol))):
    
        index_atoms_in_bonds_to_remove=locate_bonds_of_interest(mol)[i][1]
    
        index_atoms_in_cyclohxn_double_bond=locate_bonds_of_interest(mol)[i][2]
        
        bonds_to_change=locate_bonds_of_interest(mol)[i][0]
        
        bonds_to_change.append(mol.GetBondBetweenAtoms(index_atoms_in_cyclohxn_double_bond[0],index_atoms_in_cyclohxn_double_bond[1]))
                   
        copy_mol=RWMol(mol)           
            
        for bond in bonds_to_change:
                
            if bond.GetIdx()!= bonds_to_change[0].GetIdx():
                    
                copy_mol.ReplaceBond(bond.GetIdx(),double_bond)
                
            else:
                    
                copy_mol.ReplaceBond(bond.GetIdx(),single_bond)
        
            
        
        copy_mol.RemoveBond(index_atoms_in_bonds_to_remove[0],index_atoms_in_bonds_to_remove[1])
        
        copy_mol.RemoveBond(index_atoms_in_bonds_to_remove[2],index_atoms_in_bonds_to_remove[3])
                
            
        retro_DA_adduct=RWMol.GetMol(copy_mol)
            
            
        list_retroDA_adducts.append(retro_DA_adduct)
        
    return list_retroDA_adducts 


    
def get_DA_reagents(mol):
    
    list_adducts_fragments=[]
    
    list_adducts=process_retro_Diels_Alder(mol)
    
    for i in range(len(list_adducts)):
        
        list_adducts_fragments.append([])
        
        fragments=Chem.MolToSmiles(list_adducts[i]).split('.')
        
        for j in range(len(fragments)):
    
            list_adducts_fragments[i].append(Chem.MolFromSmiles(fragments[j]))       
        
    return list_adducts_fragments


def get_rDA_adducts(mol):
    
    list_adducts_fragments=[]
    
    list_adducts=process_retro_Diels_Alder(mol)
    
    for i in range(len(list_adducts)):
        
        list_adducts_fragments.append([])
        
        fragments=Chem.MolToSmiles(list_adducts[i]).split('.')
        
        for j in range(len(fragments)):
    
            list_adducts_fragments[i].append(Chem.MolFromSmiles(fragments[j]))       
        
    return list_adducts_fragments            
        

def guess_initial_transition_states(mol):  # not yet completed
    
        
    ts_geom=[]
    
    template=Chem.MolFromSmiles('C=C-C')
            
    double_bond=template.GetBonds()[0]
        
    single_bond=template.GetBonds()[1]
    
    for scheme in locate_bonds_of_interest(mol):   
       
        bond1=scheme[0][0]
        bond2=scheme[0][2]
        bond6=scheme[0][1]
        bond5=mol.GetBondBetweenAtoms(scheme[1][0],scheme[1][1])
        bond3=mol.GetBondBetweenAtoms(scheme[1][2],scheme[1][3])
        bond4=mol.GetBondBetweenAtoms(scheme[2][0],scheme[2][1])
        atom1,atom2=bond1.GetBeginAtom(),bond1.GetEndAtom()
        atom3,atom4=bond2.GetEndAtom(),bond3.GetEndAtom()
        atom5,atom6=bond4.GetEndAtom(),bond5.GetEndAtom()
        
        mol_copy=RWMol(mol)
            
        for bond in [bond1,bond2,bond4,bond6]:
        
            if bond==bond1:
        
                mol_copy.ReplaceBond(bond.GetIdx(),single_bond)
            
        
            else:
            
                mol_copy.ReplaceBond(bond.GetIdx(),double_bond)
                
        ts_geom.append(mol_copy)
        
    return ts_geom
            
        
def check_intramolecular(mol):  # Verifying if the initial molecular can result from an intramolecular Diels Alder reaction
    
    dot='.'
    
    list_check=[]
    
    for adduct in process_retro_Diels_Alder(mol):
        
        count=1
    
        for item in Chem.MolToSmiles(adduct):
            
            if item==dot:
                
                count+=1
            else:
                
                count+=0
                    
        if count==2:
                
            result=False
                
        else:
                
            result=True
            
        list_check.append(result)
        
        
    return list_check 

def compute_number_of_fragments(mol):
    
    dot='.'
    
    count=1
    
    for item in Chem.MolToSmiles(mol):
            
            if item==dot:
                
                count+=1
            else:
                
                count+=0
                    
    
    return count

def get_branching_atoms_indexes(mol,substr):
    
    branch_atom_idx=[]
         
    for atom_index in substr:
    
        atom=mol.GetAtomWithIdx(atom_index)
                
        Neighbors=[x.GetIdx() for x in atom.GetNeighbors()]            
            
            
        for atom_index_neighb in Neighbors:
                
            if atom_index_neighb not in substr:
                    
                branch_atom_idx.append(atom_index_neighb)                    
                  
    return branch_atom_idx       
        
    
def count_substituents(mol):
       
    list_number_of_substituents=[]
    
    cyclo=Chem.MolFromSmarts('C1=CCCCC1')
    
    mol=Chem.RemoveHs(mol)
    
    substructures=mol.GetSubstructMatches(cyclo,useChirality=True)
        
    for substructure in substructures:
                
        list_number_of_substituents.append(len(get_branching_atoms_indexes(mol,substructure)))
        
    
    return list_number_of_substituents


def check_ring_subst(mol):
    
    cyclo=Chem.MolFromSmarts('C1=CCCCC1')
    
    substructures=mol.GetSubstructMatches(cyclo,useChirality=True)
    
    subst_type=[]
    
    for substructure in substructures:
        
        if [mol.GetAtomWithIdx(x).IsInRing()==True for x in get_branching_atoms_indexes(mol,substructure)]:
            
            subst_type.append(True)
        else:
            subst_type.append(False)
            
            
    if True in subst_type:
        
        result=1
        
    else:
        
        result=0
            
    return result

def check_MW_category(mol):
    
    MW=rdMolDescriptors.CalcExactMolWt(mol)
    
    if MW<200:
        
        result='L'
      
    elif 200<=MW<400:
        result='M'
        
    elif 400<=MW<600:
        result='H'
        
    else:
        result='VH'
        
    return result
        
def check_NCR_category(mol):
    
    cyclo=Chem.MolFromSmarts('C1=CCCCC1')
    
    NCR=len(mol.GetSubstructMatches(cyclo,useChirality=True))
    
    if NCR==1:
        
        result=1
        
    elif NCR==2:
        
        result=2
        
    elif NCR==3:
        
        result=3
        
    else:
        
        result=4
        
    return result


def check_NS_category(mol):
    
    NS=max(count_substituents(mol))
          
    return NS


def check_DA_type(mol):
       
    if (False not in check_intramolecular(mol)):
        
        result=0
    
    elif True not in check_intramolecular(mol):
        
        result=1
        
    elif (True in check_intramolecular(mol)) and (False in check_intramolecular(mol)):
        
        result=2
        
    else:
        
        None
        
    return result
           
            
            
    
    
            
        
    
    