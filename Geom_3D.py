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


# -*- coding: utf-8 -*-


from rdkit import Chem
from Retro_DA import locate_bonds_of_interest,process_retro_Diels_Alder,check_NCR_category,get_rDA_adducts
from rdkit.Chem import Draw,AllChem,rdchem,RWMol,rdDistGeom
from rdkit.Chem.Draw import DrawingOptions
import numpy as np
import os
import shutil
from pathlib import Path
import random


"""

                4  .  .  5
             .              .
          .                    .
        .                         .
      
        3                         0
        .                       . .
           .                  . .
              .             . .       
                2  .  .  1   
                
The identification of these atoms is based on the method 'locate_bonds_of_interest()'
in the Retro_HomoDA module.


"""



DrawingOptions.includeAtomNumbers=True


def Optimize_struct(mol):

    try:
        mol=Chem.AddHs(mol)
        embed=AllChem.EmbedMolecule(mol,maxAttempts=200)
        Error=False
        if embed==-1:
            embed=AllChem.EmbedMolecule(mol,useRandomCoords=True,maxAttempts=2000,enforceChirality=False)
        if embed==-1:
            confs=rdDistGeom.EmbedMultipleConfs(mol,numConfs=20,maxAttempts=2000,enforceChirality=False,useRandomCoords=True)
            if(len(confs)<1):
                confs=rdDistGeom.EmbedMultipleConfs(mol,numConfs=40,maxAttempts=2000,enforceChirality=False,useRandomCoords=True)
            if(len(confs)<1):
                confs=rdDistGeom.EmbedMultipleConfs(mol,numConfs=60,maxAttempts=2000,enforceChirality=False,useRandomCoords=True)
            if(len(confs)>=1):
                AllChem.UFFOptimizeMolecule(mol,confId=confs[random.randint(0,len(confs)-1)],maxIters=500)
            else:
                Error=True
        else:
            AllChem.UFFOptimizeMolecule(mol,maxIters=500)

    except BaseException as e:
        print('Failed')
        print(e)
    if Error==True:
        print('Failed')
        print('Unable to embed molecule')
        
    return mol

 
def Get_Atoms_of_Interest(mol):
    
    system=mol

    paths=locate_bonds_of_interest(system)
    
    paths_atoms_of_interest=[]
    
    for path in paths:  
    
        double_bond=path[0][0]

        set_of_four_atoms_indices=path[1]

        atom_0=double_bond.GetBeginAtom()

        atom_1=double_bond.GetEndAtom()

        atom_2=system.GetAtomWithIdx(set_of_four_atoms_indices[2])

        atom_3=system.GetAtomWithIdx(set_of_four_atoms_indices[3])

        atom_4=system.GetAtomWithIdx(set_of_four_atoms_indices[1])

        atom_5=system.GetAtomWithIdx(set_of_four_atoms_indices[0])
        
        
        paths_atoms_of_interest.append([atom_0,atom_1,atom_2,atom_3,atom_4,atom_5])

    return paths_atoms_of_interest



def Gen_gjf_file_ts(mol,path,index,level,D=0.21):
    
    
    try:
        
        if os.getcwd()!=path:
        
            os.chdir(path)
            
        system=Optimize_struct(mol)
                    
                       
        for i in range(len(Get_Atoms_of_Interest(mol))):

            path=Get_Atoms_of_Interest(mol)[i] 
            
            path_number=i+1    
            
            filename='pgtsRx'+str(index)+'p'+str(path_number)
                    
            with open(filename+'.com','w') as file:
                file.write('%nprocshared=4\n')
                file.write('# opt=(modredundant,loose) '+level+' scf=xqc \n\n')
                file.write('Opt\n\n')
                file.write('0 1\n')              
    
                atom_2_index=path[2].GetIdx()
                atom_3_index=path[3].GetIdx()
                atom_4_index=path[4].GetIdx()
                atom_5_index=path[5].GetIdx()   
                    
        
                atoms=system.GetAtoms()
                        
                for atom in atoms:
                    
                    x,y,z=np.array(system.GetConformer().GetAtomPosition(atom.GetIdx()))
                    
                    x=round(x,9)
                    
                    y=round(y,9)
                    
                    z=round(z,9)
                
                    file.writelines([str(atom.GetSymbol()),'\t\t',str(x).rjust(15),'\t\t',str(y).rjust(15),'\t\t',str(z).rjust(15),'\n'])
                
                                
                file.writelines(['\n',str(atom_2_index+1),' ',str(atom_3_index+1),' S ',str(3),' ',str(D),'\n'])
                
                file.writelines([str(atom_4_index+1),' ',str(atom_5_index+1),' S ',str(3),' ',str(D),'\n']) 
                    
        
                
    except BaseException as e:
                
        print(e)               
     

    
def Gen_gjf_file_reagents(mol,path,index,level):
    
    try:
        
    
        if os.getcwd()!=path:
        
            os.chdir(path)
            
        for i in range(len(get_rDA_adducts(mol))):
                
            path_number=i+1
                
            reagents=get_rDA_adducts(mol)[i]
                
            for j in range(len(reagents)):
                    
                current_reagent=reagents[j]
                    
                filename='Rx'+str(index)+'p'+str(path_number)+'r'+str(j+1)
                    
                with open(filename+'.com','w') as file:
                        
                    file.write('%nprocshared=4\n')
                        
                    file.write('# opt freq '+level+' scf=xqc \n\n')
            
                    file.write('Optimization\n\n')
            
                    file.write('0 1\n')
            
                    Opt_mol=Optimize_struct(current_reagent)
            
                    atoms=Opt_mol.GetAtoms()
                            
                    for atom in atoms:
                        
                        x,y,z=np.array(Opt_mol.GetConformer().GetAtomPosition(atom.GetIdx()))
                        
                        x=round(x,9)
                        
                        y=round(y,9)
                        
                        z=round(z,9)
                    
                        file.writelines([str(atom.GetSymbol()),'\t\t',str(x).rjust(15),'\t\t',str(y).rjust(15),'\t\t',str(z).rjust(15),'\n'])

                    file.write('\n')
            
       
        
    except BaseException as e:
        
        print(e)
    
    
    
def Gen_gjf_file_adducts(mol,path,index,level):    
    
    try:
        
        os.chdir(path)
        
        filename='Cycloadduct'+str(index)
        
        with open(filename+'.com','w') as file:
            
            file.write('%nprocshared=4\n')
            
            file.write('# opt freq '+level+' scf=xqc \n\n')
            
            file.write('Optimization\n\n')
            
            file.write('0 1\n')
            
            Opt_mol=Optimize_struct(mol)
            
            atoms=Opt_mol.GetAtoms()
                            
            for atom in atoms:
                        
                x,y,z=np.array(Opt_mol.GetConformer().GetAtomPosition(atom.GetIdx()))
                        
                x=round(x,9)
                        
                y=round(y,9)
                        
                z=round(z,9)
                    
                file.writelines([str(atom.GetSymbol()),'\t\t',str(x).rjust(15),'\t\t',str(y).rjust(15),'\t\t',str(z).rjust(15),'\n'])
                
            file.write('\n')
             
    except BaseException as e:
        
        print(e)
                                
    
def Gen_AOI_file(cwd):
    
    import sys    
    
    
    if os.path.isfile(os.path.join(cwd,'SMILES.txt')):
        smiles=open(os.path.join(cwd,'SMILES.txt'),'r')
        content=smiles.readlines()
        smiles.close()
    else:
        sys.exit('SMILES.txt file not found !!! Cannot AOI.csv file')
    
    with open('AOI.csv','w') as file:
    
        file.write('Reaction_ID,Atom0,Atom1,Atom2,Atom3,Atom4,Atom5\n')
        
        for r in range(len(content)):
        
            mol=Chem.MolFromSmiles(content[r])

            for i in range(len(Get_Atoms_of_Interest(mol))):

                atoms_cp=Get_Atoms_of_Interest(mol)[i] 

                path_number=i+1

                ID='Rx'+str(r+1)+'p'+str(path_number)

                atoms=[str(x.GetIdx()+1) for x in atoms_cp]

                file.writelines([ID,',',atoms[0],',',atoms[1],',',atoms[2],',',atoms[3],',',atoms[4],',',atoms[5],'\n'])



