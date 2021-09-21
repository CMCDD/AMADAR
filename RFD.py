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
from IRC import IRC_1F,IRC_2F
from TS import Opt_TS



class RF_atomic_decomp():
    
    def __init__(self,filename):
        
        self.filename=filename

    def get_list_atomic_symbols(self):
        
        geom=IRC_1F(self.filename).Df_geoms_string()[0]
        
        return geom['Atomic Symbol'].tolist()

    def get_atomic_force_string(self,atomic_index):
        
        filename=self.filename
    
        force_matrices=IRC_1F(filename).Df_force_matrix_string()
        
        list_af=[]
        
        count=0
        
        for matrix in force_matrices:
            
            list_af.append([])
            
            list_af[count].append(matrix['Fx'][atomic_index-1])
            list_af[count].append(matrix['Fy'][atomic_index-1])
            list_af[count].append(matrix['Fz'][atomic_index-1])
            
            count+=1
            
        return list_af
    
    def get_atomic_position_string(self,atomic_index):
        
        filename=self.filename
        
        geoms=IRC_1F(filename).Df_geoms_string()
        
        list_atm_pos=[]
        
        count=0
        
        for geom in geoms:
            
            list_atm_pos.append([])
            
            list_atm_pos[count].append(geom['X'][atomic_index-1])
            list_atm_pos[count].append(geom['Y'][atomic_index-1])
            list_atm_pos[count].append(geom['Z'][atomic_index-1])
            
            count+=1
            
        return list_atm_pos
    
    def get_atomic_pos_vect_string(self,atomic_index):
        
        import numpy as np
        
        positions=self.get_atomic_position_string(atomic_index)
        
        origin=positions[0]
        
        list_atom_vect=[]
        
        for position in positions:
            
            actual_pos_vec=(np.array(position)-np.array(origin)).tolist()
            
            list_atom_vect.append(actual_pos_vec)
            
            
        return list_atom_vect
    
    
    def derivation(self,list_y,list_x):
        
        
        first_term=(list_y[1]-list_y[0])/(list_x[1]-list_x[0])
        
        last_term=(list_y[-2]-list_y[-1])/(list_x[-2]-list_x[-1])
        
        list_der=[first_term]
        
        if len(list_y)==len(list_x):
            
            for j in range(1,len(list_y)-1):
                
                back_der=(list_y[j+1]-list_y[j])/(list_x[j+1]-list_x[j])
                
                forw_der=(list_y[j]-list_y[j-1])/(list_x[j]-list_x[j-1])
                
                der=(back_der+forw_der)*0.5
                
                list_der.append(der)
                
        else:
            
            print('List lengths do not match during derivation attempt')
            
        list_der.append(last_term)
        
        
        return list_der
    
    def pos_vect_grad_string(self,atomic_index):
        
        
        pos_vect_string=self.get_atomic_pos_vect_string(atomic_index)
        
        Rx_coords=IRC_1F(self.filename).get_Rx_coords()[0]
        
        pos_vect_string_x_comp=[r[0] for r in pos_vect_string ]
        
        pos_vect_string_y_comp=[r[1] for r in pos_vect_string ]
        
        pos_vect_string_z_comp=[r[2] for r in pos_vect_string ]
        
        
        der_x=self.derivation(pos_vect_string_x_comp,Rx_coords)
        der_y=self.derivation(pos_vect_string_y_comp,Rx_coords)
        der_z=self.derivation(pos_vect_string_z_comp,Rx_coords)
        
        der=[(x,y,z) for x,y,z in zip(der_x,der_y,der_z)]
        
        return der, der_x, der_y, der_z
    
    
    def get_atomic_reaction_force(self,atomic_index):
        
        import numpy as np
        
        #import matplotlib.pyplot as plt
        
        atomic_forces=self.get_atomic_force_string(atomic_index)
        
        pos_vect_grad=self.pos_vect_grad_string(atomic_index)[0]
        
        atom_rf=[]
        
        Rx_coords=IRC_1F(self.filename).get_Rx_coords()[0]
        
        if len(atomic_forces)==len(pos_vect_grad):
            
            for i in range(len(atomic_forces)):
                
                atom_rf.append(sum((np.array(atomic_forces[i])*np.array(pos_vect_grad[i])).tolist()))
            
        else:
            
            print('Vector lenghts do not match while computing atomic reaction force')
        
        
        #plt.plot(Rx_coords,atom_rf)
        #plt.show()
        
        return atom_rf
    
    def get_fragment_reaction_force(self,list_atomic_index_frag):
        
        frag_arf=[]
        
        atom_rfs=[]
        
        
        for atomic_index in list_atomic_index_frag:
            
                        
            atom_rf=self.get_atomic_reaction_force(atomic_index)
            
            atom_rfs.append(atom_rf)
            
            
        for i in range(len(atom_rfs[0])):
            
            sum_rfs=0
            
            for j in range(len(atom_rfs)):
                
                sum_rfs+=atom_rfs[j][i]
                
            
            frag_arf.append(sum_rfs)
        
            
        
        return frag_arf               
               
    
    def plot_multiple_atomic_decomp(self,list_atomic_index,output_dir):
        
        import matplotlib.pyplot as plt
        
        Rx_coords=IRC_1F(self.filename).get_Rx_coords()[0]
        
        alpha_limit=IRC_1F(self.filename).get_min_max_Rx_force()[0]
        gamma_limit=IRC_1F(self.filename).get_min_max_Rx_force()[1]
        atomic_symbols=self.get_list_atomic_symbols()
        
        
            
        plt.figure(figsize=(6,4),dpi=500)
        
        plt.axvline(alpha_limit,c='black',linestyle='--')
        plt.axvline(gamma_limit,c='black',linestyle='--')
        plt.axvline(0,c='pink',alpha=0.4)
        plt.axhline(0,c='black',linestyle='--',alpha=0.4)
        plt.axvspan(min(Rx_coords),alpha_limit,alpha=0.3,color='steelblue')        
        plt.axvspan(gamma_limit,max(Rx_coords),alpha=0.3,color='steelblue')
        plt.ylabel('Force (hartree/bohr)',labelpad=8)
        plt.xlabel('Reaction coordinate (ξ)',labelpad=8)
        
        #Creating symmetrical plot boundaries
        if abs(Rx_coords[-1])<abs(Rx_coords[0]):
            plt.xlim(-Rx_coords[-1],Rx_coords[-1])
        else:
            plt.xlim(Rx_coords[0],-Rx_coords[0])
            
        #Plotting atomic components
        
        for atomic_index in list_atomic_index:
            
            atom_rf=self.get_atomic_reaction_force(atomic_index)
            
            plt.plot(Rx_coords,atom_rf,label=atomic_symbols[atomic_index-1]+str(atomic_index))
            
        plt.legend(loc='best',fontsize=9)
        
        tail=os.path.split(self.filename)[1]
        figname=os.path.join(output_dir,'rfad'+tail[3:-3]+'png')
        plt.savefig(figname)
        
    def get_rf_from_atom_contr(self):
        
        #Computing the total reaction force from atomic contributions
        
        mol_size=len(IRC_1F(self.filename).Df_geoms_string()[0])
        
        all_atom_ind=[(x+1) for x in range(mol_size)]
        
        rf=self.get_fragment_reaction_force(all_atom_ind)
        
        return rf
        
    def plot_multiple_fragments_force_decomp(self,list_fragments,output_dir,list_fragment_names=None):
        
        import matplotlib.pyplot as plt
        
        #Default group names
        
        default=["Group "+str(x+1) for x in range(len(list_fragments))]
        
        if list_fragment_names==None:
            
            list_fragment_names=default
        
        else:
            
            list_fragment_names=list_fragment_names
            
        #Computing fragment contributions to the reaction force
        
        frag_arfs=[]
        
        for frag in list_fragments:
            
            arf=self.get_fragment_reaction_force(frag)
            
            frag_arfs.append(arf)
            
        #Plotting
        
        Rx_coords=IRC_1F(self.filename).get_Rx_coords()[0]        
        alpha_limit=IRC_1F(self.filename).get_min_max_Rx_force()[0]
        gamma_limit=IRC_1F(self.filename).get_min_max_Rx_force()[1]  
        reaction_force=self.get_rf_from_atom_contr()
        
        plt.figure(figsize=(6,4),dpi=500)        
        plt.axvline(alpha_limit,c='black',linestyle='--')
        plt.axvline(gamma_limit,c='black',linestyle='--')
        plt.axvline(0,c='pink',alpha=0.4)
        plt.axhline(0,c='black',linestyle='--',alpha=0.4)
        plt.axvspan(min(Rx_coords),alpha_limit,alpha=0.3,color='steelblue')        
        plt.axvspan(gamma_limit,max(Rx_coords),alpha=0.3,color='steelblue')
        plt.ylabel('Force (hartree/bohr)',labelpad=8)
        plt.xlabel('Reaction coordinate (ξ)',labelpad=8)
        plt.plot(Rx_coords,reaction_force,linestyle='dotted',color='black',label='F')
        
        #Creating symmetrical plot boundaries
        if abs(Rx_coords[-1])<abs(Rx_coords[0]):
            plt.xlim(-Rx_coords[-1],Rx_coords[-1])
        else:
            plt.xlim(Rx_coords[0],-Rx_coords[0])
            
                
        #Plotting fragment components
        i=0
        for frag_rf in frag_arfs:
            
            plt.plot(Rx_coords,frag_rf,label='F('+list_fragment_names[i]+')')
            i+=1
            
        plt.legend(loc='best',fontsize=8)
        
        tail=os.path.split(self.filename)[1]
        figname=os.path.join(output_dir,'rffd'+tail[3:-3]+'png')
        plt.savefig(figname)
        
        
    
            
class RFC_atomic_decomp(RF_atomic_decomp):

    def __init__(self,filename):
        
        super().__init__(filename)   #inheriting methods for the RF_atomic_decomp() class
        
        self.filename=filename        
        
        
    def get_atomic_reaction_force_constant(self,atomic_index):
        
        rf=self.get_atomic_reaction_force(atomic_index)
        
        Rx_coords=IRC_1F(self.filename).get_Rx_coords()[0]
        
        rfc=[-x for x in self.derivation(rf,Rx_coords)]
        
        return rfc
    
    def get_rfc_from_atom_contr(self):
        
        #Computing the total reaction force constant from atomic contributions
        
        mol_size=len(IRC_1F(self.filename).Df_geoms_string()[0])
        
        all_atom_ind=[(x+1) for x in range(mol_size)]
        
        rfc=self.get_fragment_reaction_force_constant(all_atom_ind)
        
        return rfc
    
    def plot_multiple_atomic_decomp(self,list_atomic_index,output_dir):
        
        import matplotlib.pyplot as plt
        
        Rx_coords=IRC_1F(self.filename).get_Rx_coords()[0]
        
        alpha_limit=IRC_1F(self.filename).get_min_max_Rx_force()[0]
        gamma_limit=IRC_1F(self.filename).get_min_max_Rx_force()[1]
        atomic_symbols=self.get_list_atomic_symbols()
        
        
            
        plt.figure(figsize=(6,4),dpi=500)
        
        plt.axvline(alpha_limit,c='black',linestyle='--')
        plt.axvline(gamma_limit,c='black',linestyle='--')
        plt.axvline(0,c='pink',alpha=0.4)
        plt.axhline(0,c='black',linestyle='--',alpha=0.4)
        plt.axvspan(min(Rx_coords),alpha_limit,alpha=0.3,color='steelblue')        
        plt.axvspan(gamma_limit,max(Rx_coords),alpha=0.3,color='steelblue')
        plt.ylabel('Force constant (hartree/$bohr^2$)',labelpad=8)
        plt.xlabel('Reaction coordinate (ξ)',labelpad=8)
        
        #Creating symmetrical plot boundaries
        if abs(Rx_coords[-1])<abs(Rx_coords[0]):
            plt.xlim(-Rx_coords[-1],Rx_coords[-1])
        else:
            plt.xlim(Rx_coords[0],-Rx_coords[0])
            
        #Plotting atomic components
        
        for atomic_index in list_atomic_index:
            
            atom_rfc=self.get_atomic_reaction_force_constant(atomic_index)
            
            plt.plot(Rx_coords,atom_rfc,label=atomic_symbols[atomic_index-1]+str(atomic_index))
            
        plt.legend(loc='best',fontsize=9)
        
        tail=os.path.split(self.filename)[1]
        figname=os.path.join(output_dir,'rfcad'+tail[3:-3]+'png')
        plt.savefig(figname)
        
        
    def get_fragment_reaction_force_constant(self,list_atomic_index_frag):
        
        frag_rf=self.get_fragment_reaction_force(list_atomic_index_frag)
        
        Rx_coords=IRC_1F(self.filename).get_Rx_coords()[0]
        
        frag_rfc=[-x for x in self.derivation(frag_rf,Rx_coords)]
            
        
        return frag_rfc
    
    
    def plot_multiple_fragments_force_constant_decomp(self,list_fragments,output_dir,list_fragment_names=None):
        
        import matplotlib.pyplot as plt
        
        #Default group names
        
        default=["Group "+str(x+1) for x in range(len(list_fragments))]
        
        if list_fragment_names==None:
            
            list_fragment_names=default
        
        else:
            
            list_fragment_names=list_fragment_names
            
        #Computing fragment contributions to the reaction force
        
        frag_arfcs=[]
        
        for frag in list_fragments:
            
            arfc=self.get_fragment_reaction_force_constant(frag)
            
            frag_arfcs.append(arfc)
            
        #Plotting
        
        Rx_coords=IRC_1F(self.filename).get_Rx_coords()[0]        
        alpha_limit=IRC_1F(self.filename).get_min_max_Rx_force()[0]
        gamma_limit=IRC_1F(self.filename).get_min_max_Rx_force()[1]  
        reaction_force_constant=self.get_rfc_from_atom_contr()
        
        plt.figure(figsize=(6,4),dpi=500)        
        plt.axvline(alpha_limit,c='black',linestyle='--')
        plt.axvline(gamma_limit,c='black',linestyle='--')
        plt.axvline(0,c='pink',alpha=0.4)
        plt.axhline(0,c='black',linestyle='--',alpha=0.4)
        plt.axvspan(min(Rx_coords),alpha_limit,alpha=0.3,color='steelblue')        
        plt.axvspan(gamma_limit,max(Rx_coords),alpha=0.3,color='steelblue')
        plt.ylabel('Force constant (hartree/$bohr^2$)',fontsize=10,labelpad=8)
        plt.xlabel('Reaction coordinate (ξ)',fontsize=10,labelpad=8)
        plt.plot(Rx_coords,reaction_force_constant,linewidth=2,linestyle='dashed',color='black',label='κ(atoms)')
        
        #Creating symmetrical plot boundaries
        if abs(Rx_coords[-1])<abs(Rx_coords[0]):
            plt.xlim(-Rx_coords[-1],Rx_coords[-1])
        else:
            plt.xlim(Rx_coords[0],-Rx_coords[0])
            
                
        #Plotting fragment components
        i=0
        colors=['violet','maroon']
        for frag_rfc in frag_arfcs:
            
            plt.plot(Rx_coords,frag_rfc,c=colors[i],label='κ('+list_fragment_names[i]+')')
            i+=1
            
        plt.legend(loc='best',fontsize=8)
        tail=os.path.split(self.filename)[1]
        figname=os.path.join(output_dir,'rfcfd'+tail[3:-3]+'png')
        plt.savefig(figname)   
            
    def plot_K_ATOMS_BONDS(self):
            
        import matplotlib.pyplot as plt
            
        rfc_atoms=self.get_rfc_from_atom_contr()
            
        rfc_all=[round(x/627.509,4) for x in IRC_1F(self.filename).get_Rx_force_constant()[0]]
            
        rfc_bonds=[(x-y)for x,y in zip(rfc_all,rfc_atoms)]
        
        Rx_coords=IRC_1F(self.filename).get_Rx_coords()[0]        
        alpha_limit=IRC_1F(self.filename).get_min_max_Rx_force()[0]
        gamma_limit=IRC_1F(self.filename).get_min_max_Rx_force()[1]  
        reaction_force_constant=self.get_rfc_from_atom_contr()
        
        plt.figure(figsize=(6,4),dpi=1500)        
        plt.axvline(alpha_limit,c='black',linestyle='--')
        plt.axvline(gamma_limit,c='black',linestyle='--')
        plt.axvline(0,c='pink',alpha=0.4)
        plt.axhline(0,c='black',linestyle='--',alpha=0.4)
        plt.axvspan(min(Rx_coords),alpha_limit,alpha=0.3,color='teal')        
        plt.axvspan(gamma_limit,max(Rx_coords),alpha=0.3,color='teal')
        plt.ylabel('Force constant (hartree/$bohr^2$)',labelpad=8)
        plt.xlabel('Reaction coordinate (ξ)',labelpad=8)
        
        #Creating symmetrical plot boundaries
        if abs(Rx_coords[-1])<abs(Rx_coords[0]):
            plt.xlim(-Rx_coords[-1],Rx_coords[-1])
        else:
            plt.xlim(Rx_coords[0],-Rx_coords[0])
            
        #Plotting   
            
        plt.plot(Rx_coords,rfc_atoms,linestyle='dotted',color='tomato',label='κ(atoms)')
        plt.plot(Rx_coords,rfc_bonds,linestyle='dashed',color='violet',label='κ(bonds)')
        plt.plot(Rx_coords,rfc_all,linestyle='dotted',color='black',label='κ')
        plt.legend(loc='best',fontsize=8)
            
            