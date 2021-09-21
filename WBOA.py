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

import os
import shutil
from tqdm import tqdm
from TS import Opt_TS
from IRC import IRC_1F


dict_atoms={
    
    1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C',
    
    7: 'N', 8:'O', 9:'F', 10:'Ne', 11:'Na', 12: 'Mg',
    
    13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl', 18: 'Ar',
    
    19: 'K', 20: 'Ca', 21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',
    
    27:'Co',28:'Ni',29:'Cu',30:'Zn',31:'Ga',32:'Ge',33:'As',34:'Se',
    
    35:'Br',36:'Kr',37:'Rb',38:'Sr',39:'Y',40:'Zr',41:'Nb',42:'Mo',43:'Tc',
    
    44:'Ru',45:'Rh',46:'Pd',47:'Ag',48:'Cd',49:'In',50:'Sn',51:'Sb',52:'Te',
    
    53:'I',54:'Xe'}


class WBO():
    
    def __init__(self,file):
        
        self.filename=file 
        
        
    def get_AOI_table(self):
        
        import pandas as pd
        
        file=self.filename
        
        file_tail=os.path.split(file)[1]
        
        AOI_table=pd.read_csv('AOI.csv')
        
        rx_IDs=AOI_table['Reaction_ID']
        
        rx_IDs=rx_IDs.tolist()
        
        return AOI_table,rx_IDs
    
    
    def get_atoms_of_interest(self):
        
        file=self.filename
        
        file_tail=os.path.split(file)[1]
        
        AOI_table=self.get_AOI_table()[0]
        
        rx_IDs=self.get_AOI_table()[1]
        
        current_Rx_ID=file_tail.split('pt')[0][3:]
        
        current_Rx_index=rx_IDs.index(current_Rx_ID)
        
        atom0,atom1=AOI_table['Atom0'][current_Rx_index],AOI_table['Atom1'][current_Rx_index]
        
        atom2,atom3=AOI_table['Atom2'][current_Rx_index],AOI_table['Atom3'][current_Rx_index]
        
        atom4,atom5=AOI_table['Atom4'][current_Rx_index],AOI_table['Atom5'][current_Rx_index]
                
            
        return atom0,atom1,atom2,atom3,atom4,atom5
        
    
    def locate_WBO_matrix(self):
        
        file=self.filename
        
        ref_word1='Wiberg bond index matrix in the NAO basis:     '
        
        ref_word2='Atom'
        
        ref_word3='Wiberg bond index, Totals by atom:'
        
        myfile=open(file,'r')
        
        nl=len(myfile.readlines())
        
        myfile.seek(0,0)
        
        WBO_matrix=[]
        
        for i in range(nl):
            
            lc=myfile.readline()
            
            if ref_word1 in lc:
                
                begin=myfile.tell()-len(lc)
                
                
            if ref_word2 in lc: 
                
                WBO_matrix.append(myfile.tell()-len(lc))
                
                                
            if ref_word3 in lc:
                
    
                break
        result=[]
        
        for value in WBO_matrix:
            
                        
            if value>=begin:
                
                result.append(value)
        
        
        myfile.close()
        
        
        return result        
    
    
    def get_WBO_parts(self):
        
        import pandas as pd 
        
        
        file=self.filename
        
        pos=self.locate_WBO_matrix()
        
        
        WBO_matrix_parts=[]
        
        with open(file,'r') as myfile:
            
            for value in pos:
                
                matrix_part=[]
                
                myfile.seek(value,0)
                    
                Columns=myfile.readline().split()[1:]
                
                myfile.readline()
                
                count_lines=0
                
                while True:
                    
                    lc=myfile.readline()
                    
                    if len(lc.split())>=2:
                        
                        matrix_part.append([])
                        
                        for j in range(len(lc.split())-2):
                            
                            matrix_part[count_lines].append(lc.split()[j+2])
                            
                        count_lines+=1
                        
                    else:
                        
                        break
                    
                data=pd.DataFrame(matrix_part,columns=Columns)
                
                WBO_matrix_parts.append(data)
                
                
        return WBO_matrix_parts
                            

    def build_WBO_matrix(self):
        
        import pandas as pd       
        
        
        WBO_matrix_parts=self.get_WBO_parts()
        
        
        WBO_matrix=pd.DataFrame()
        
        WBO_matrix['Atoms']=[x for x in range(1,len(WBO_matrix_parts[0])+1)]        
        
        for part in WBO_matrix_parts:
            
            WBO_matrix=pd.concat([WBO_matrix,part],axis=1)
            
            
        return WBO_matrix
   
    def get_bonding_atoms(self):      
        
        AOI=self.get_atoms_of_interest()
        atom0,atom1,atom2=AOI[0],AOI[1],AOI[2]
        atom3,atom4,atom5=AOI[3],AOI[4],AOI[5]
        BA=['C'+str(atom5)+'-'+'C'+str(atom0),'C'+str(atom0)+'-'+'C'+str(atom1),'C'+str(atom1)+'-'+'C'+str(atom2)
            ,'C'+str(atom2)+'-'+'C'+str(atom3),'C'+str(atom3)+'-'+'C'+str(atom4),'C'+str(atom4)+'-'+'C'+str(atom5)]
        
        return BA
    
    def get_bond_orders(self):
        
        AOI=self.get_atoms_of_interest()
        atom0,atom1,atom2=AOI[0],AOI[1],AOI[2]
        atom3,atom4,atom5=AOI[3],AOI[4],AOI[5]
        WBO_matrix=self.build_WBO_matrix()
        
        BO_5_0=WBO_matrix.loc[atom5-1,str(atom0)]
        BO_0_1=WBO_matrix.loc[atom0-1,str(atom1)]
        BO_1_2=WBO_matrix.loc[atom1-1,str(atom2)]
        BO_2_3=WBO_matrix.loc[atom2-1,str(atom3)]
        BO_3_4=WBO_matrix.loc[atom3-1,str(atom4)]
        BO_4_5=WBO_matrix.loc[atom4-1,str(atom5)]
        
        bonds=[BO_5_0,BO_0_1,BO_1_2,BO_2_3,BO_3_4,BO_4_5]
        
        
        return bonds
    
class WBO_analysis_1F():
    
    def __init__(self,list_files,irc_file):
        
        self.files=list_files
        
        self.irc=irc_file
        
        
    
    def sort_files(self):
        
        files=self.files 
        
        nf=len(files)
        
        sorted_list=[]
        
        for i in range(1,nf+1):
            
            for file in files:
                
                tail=os.path.split(file)[1]
                
                if int(tail.split('pt')[1][:-4])==i:
                    
                    sorted_list.append(file)
                    
                    break
                
                
        return sorted_list
    
    def reverse_list(self,my_list):
        
        new_list=[]
        
        for i in range(len(my_list)):
            
            new_list.append(my_list[len(my_list)-1-i])
            
        return new_list
    
    def divide_IRC(self):
        
        irc=IRC_1F(self.irc)
        
        alpha=irc.get_min_max_Rx_force()[0]
        
        gamma=irc.get_min_max_Rx_force()[1]
        
        return alpha,gamma
    
    def plot(self,output_dir):
        
        import matplotlib.pyplot as plt
        
        alpha_Rc,gamma_Rc=self.divide_IRC()[0],self.divide_IRC()[1]
        
        files=self.sort_files()             #files=self.reverse_list(files)
        
        output_filename=os.path.join(output_dir,os.path.split(self.irc)[1])
        
        BOs=[]        
        
        for file in files:
            
            job=WBO(file)
            
            curr_WBOs=job.get_bond_orders()
            
            BOs.append(curr_WBOs)
            
        
        bond_5_0_string=[eval(x[0]) for x in BOs]
        bond_0_1_string=[eval(x[1]) for x in BOs]
        bond_1_2_string=[eval(x[2]) for x in BOs]
        bond_2_3_string=[eval(x[3]) for x in BOs]
        bond_3_4_string=[eval(x[4]) for x in BOs]
        bond_4_5_string=[eval(x[5]) for x in BOs]
        
        
        BA=WBO(files[0]).get_bonding_atoms()
        
        Rx_coords=IRC_1F(self.irc).get_Rx_coords()[0]
        RFC=IRC_1F(self.irc).get_Rx_force_constant()[0]
        RFC_max_index=(RFC).index(max(RFC))
        
        max_pos=Rx_coords[RFC_max_index]
        
        plt.figure(figsize=(7,5))
        #Creating symmetrical plot boundaries
        if abs(Rx_coords[-1])<abs(Rx_coords[0]):
            plt.xlim(-Rx_coords[-1],Rx_coords[-1]+2)
        else:
            plt.xlim(Rx_coords[0],-Rx_coords[0]+2)   
            
        legend_properties = {'weight':'bold'}       
        plt.plot(Rx_coords,bond_5_0_string,linewidth=2,c='blue',label=BA[0])
        plt.plot(Rx_coords,bond_0_1_string,linewidth=2,c='red',label=BA[1])
        plt.plot(Rx_coords,bond_1_2_string,linewidth=2,linestyle='--',c='green',label=BA[2])
        plt.plot(Rx_coords,bond_2_3_string,linewidth=2,c='yellow',label=BA[3])
        plt.plot(Rx_coords,bond_3_4_string,linewidth=2,c='cyan',label=BA[4])
        plt.plot(Rx_coords,bond_4_5_string,linewidth=2,linestyle='--',c='magenta',label=BA[5])
        plt.axvline(alpha_Rc,c='black',linestyle='--',alpha=0.75)
        plt.axvline(gamma_Rc,c='black',linestyle='--',alpha=0.75)
        plt.axhline(2,linestyle='dashdot',c='black')
        plt.axhline(1,linestyle='dashdot',c='black')
        plt.axvspan(Rx_coords[0]-1,alpha_Rc,color='silver',alpha=0.5)
        plt.axvspan(gamma_Rc,8,color='silver',alpha=0.5)
        plt.axvline(0,c='pink',linestyle='--')
        plt.axvline(max_pos,c='black')
        plt.ylabel('Wiberg Bond Order',fontsize=8,fontweight='bold',labelpad=8)
        plt.xlabel('Reaction coordinate (ξ)',fontsize=8,labelpad=10,fontweight='bold')
        plt.text(max_pos+0.4,1.5, 'ω',fontweight='bold',fontsize=12)
        plt.legend(loc='best',fontsize=9,prop=legend_properties)
        
        
        plt.savefig(output_filename[:-3]+'png',dpi=500)
        
        
            