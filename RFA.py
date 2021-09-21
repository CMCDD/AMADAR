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

dict_atoms={
    
    1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C',
    
    7: 'N', 8:'O', 9:'F', 10:'Ne', 11:'Na', 12: 'Mg',
    
    13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl', 18: 'Ar',
    
    19: 'K', 20: 'Ca', 21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',
    
    27:'Co',28:'Ni',29:'Cu',30:'Zn',31:'Ga',32:'Ge',33:'As',34:'Se',
    
    35:'Br',36:'Kr',37:'Rb',38:'Sr',39:'Y',40:'Zr',41:'Nb',42:'Mo',43:'Tc',
    
    44:'Ru',45:'Rh',46:'Pd',47:'Ag',48:'Cd',49:'In',50:'Sn',51:'Sb',52:'Te',
    
    53:'I',54:'Xe'}


class RFA_1file():
    
    def __init__(self,file):
        
        self.filename=file
        
    def get_number_of_lines(self):
        
        file=self.filename

        irc=open(file)

        nl=len(irc.readlines())

        irc.close()
        
        return nl
    
    def reverse_list(self,my_list):
        
        new_list=[]
        
        for i in range(len(my_list)):
            
            new_list.append(my_list[len(my_list)-1-i])
            
        return new_list
    
    def get_Rx_energies_VS_coordinates(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Energies reported relative to the TS energy'
        
        pos=None
        
        list_rx_coords_ener=[]
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                    
                        pos=my_file.tell()-len(lc)
                        
                        break
                    
                my_file.seek(pos,0)
                
                for i in range(5):
                    my_file.readline()
                
                
                count=0
                while True:
                    
                    lc=my_file.readline()
                    
                    if '-----------------' not in lc:
                        
                        list_rx_coords_ener.append([])
                        
                        list_rx_coords_ener[count].append(eval(lc.split()[2]))
                        
                        list_rx_coords_ener[count].append(eval(lc.split()[1]))
                        
                        count+=1
                        
                        #print(lc)
                        
                    else:
                        
                        break
                        
                        
            
        else:
            
            print('File Not Found')
            
            
        return list_rx_coords_ener
    
    def get_Rx_coords(self):
        
        E_Rc=self.get_Rx_energies_VS_coordinates()
        
        Rc=[]
        
        E=[]
        
        for i  in range(len(E_Rc)):
            
            Rc.append(E_Rc[i][0])
            
            E.append(E_Rc[i][1])
            
        if abs(E[0])>abs(E[-1]) and Rc[0]<Rc[-1]:
            
            Rc=[-x for x in self.reverse_list(Rc)]
            
        if abs(E[0])>abs(E[-1]) and Rc[0]>Rc[-1]:
            
            Rc=self.reverse_list(Rc)
        
        
        return Rc
    
    def get_Energy(self):
        
        E_Rx_coords=self.get_Rx_energies_VS_coordinates()
        
        E=[]
        
        for i in range(len(E_Rx_coords)):
            
            E.append(E_Rx_coords[i][1])
            
        if abs(E[-1])>abs(E[0]):
            
            E=[(x-E[0])*627.509 for x in E]
            
            order='Normal'
        else:
            
            E=[(x-E[-1])*627.509 for x in E]
            
            E=self.reverse_list(E)
            
            order='Reverse'
            
        
        return E,order        
    
    def get_min_max_Rx_force(self):
        
        force=self.get_Rx_force()[0]
        
        Rc=self.get_Rx_coords()
        
        global_min=Rc[0]
        
        global_max=Rc[0]
        
        test=force[0]
        
        for i in range(1,len(force)):
            
            if force[i]<test:
                
                global_min=Rc[i]
                
                test=force[i]
                
        
        test=force[0]
        
        for j in range(1,len(force)):
                
            if force[j]>test:
                
                global_max=Rc[j]
                
                test=force[j]
                
        
        return global_min,global_max
    
    def get_Rx_force(self):
        
        E=self.get_Energy()[0]
        
        Rc=self.get_Rx_coords()
        
        centr_der=[]
        
        first_term=(E[1]-E[0])/(Rc[1]-Rc[0])
        
        centr_der.append(first_term)
        
        for i in range(1,len(E)-1):
            
            forw_der=(E[i+1]-E[i])/(Rc[i+1]-Rc[i])
            
            backw_der=(E[i]-E[i-1])/(Rc[i]-Rc[i-1])
            
            centr_der.append((forw_der+backw_der)/2)
            
        centr_der.append((E[-2]-E[-1])/(Rc[-2]-Rc[-1]))
        
        force=[-x for x in centr_der]
        
        return force,Rc
    
    
    def get_Rx_force_constant(self):
        
        force=self.get_Rx_force()[0]
        
        Rc=self.get_Rx_coords()
        
        centr_der=[(force[1]-force[0])/(Rc[1]-Rc[0])]
        
        for i in range(1,len(force)-1):
            
            forw_der=(force[i+1]-force[i])/(Rc[i+1]-Rc[i])
            
            backw_der=(force[i]-force[i-1])/(Rc[i]-Rc[i-1])
            
            centr_der.append((forw_der+backw_der)/2)
            
        centr_der.append((force[-2]-force[-1])/(Rc[-2]-Rc[-1]))
            
        force_constant=[-x for x in centr_der]
        
        return force_constant,Rc 
    
    def get_TS_region_extend(self):
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]
        
        return round(abs(gamma_limit-alpha_limit),2)
    
    def dec_activ_energy(self):

        Alpha_index=(self.get_Rx_coords()).index(self.get_min_max_Rx_force()[0])
        Gamma_index=(self.get_Rx_coords()).index(self.get_min_max_Rx_force()[1])
        Energy=self.get_Energy()[0]
        E_TS=max(Energy)
        E_Reactant=Energy[0]
        E_Alpha=Energy[Alpha_index]
        E_Gamma=Energy[Gamma_index]
        E_Prod=Energy[-1]
        
        return round(E_Alpha-E_Reactant,2),round(E_TS-E_Alpha,2),round(E_Gamma -E_TS,2),round(E_Prod-E_Gamma,2)  
    
    def dec_react_force_energy(self):
        
        Alpha_index=(self.get_Rx_coords()).index(self.get_min_max_Rx_force()[0])
        Gamma_index=(self.get_Rx_coords()).index(self.get_min_max_Rx_force()[1])
        RF=self.get_Rx_force()[0]
        RF_TS=0
        RF_Reactant=RF[0]
        RF_Alpha=RF[Alpha_index]
        RF_Gamma=RF[Gamma_index]
        RF_Prod=RF[-1]
        
        return round(RF_Alpha-RF_Reactant,2),round(RF_TS-RF_Alpha,2),round(RF_Gamma -RF_TS,2),round(RF_Prod-RF_Gamma,2)  
        
    def plot(self,output_dir):
        
        import matplotlib.pyplot as plt
        
        
        file_tail=os.path.split(self.filename)[1]
        
        filename=os.path.join(output_dir,'rfa'+file_tail[3:-4])
        
        title=(os.path.split(self.filename)[1])[4:-6]
        
        alpha=self.get_min_max_Rx_force()[0]
        
        gamma=self.get_min_max_Rx_force()[1]
        
        alpha_index=(self.get_Rx_coords()).index(self.get_min_max_Rx_force()[0])
        
        BETA=(self.get_Energy()[0]).index(max(self.get_Energy()[0]))
        
        GAMMA=(self.get_Rx_coords()).index(self.get_min_max_Rx_force()[1])
        
        TS_extent=self.get_TS_region_extend()
        
        Energy=self.get_Energy()[0]
        
        Rx_force=self.get_Rx_force()[0]
        
        
        Rx_force_constant=self.get_Rx_force_constant()[0]
        
        Rx_coords=self.get_Rx_coords()
                
        plt.figure(figsize=(6,4))
        
        
        fig,ax=plt.subplots()
        
        ax.plot(Rx_coords,Energy,c='blue',linewidth=2,label='E')
        
        ax.set_ylabel('Energy (kcal/mol)',fontsize=12)
        
        ax.set_xlabel('Reaction coordinate (ξ)',fontsize=12,labelpad=8)
        
        ax.set_ylim(min(Energy)-5,max(Energy)+5)
        
        plt.legend(loc='upper left')
        
        plt.axhline(0,c='black',alpha=0.40)
        
        plt.axvspan(min(Rx_coords),alpha,alpha=0.3,color='silver')
        
        plt.axvspan(gamma,max(Rx_coords),alpha=0.3,color='silver')
        
        
        ax2=ax.twinx()
        
        ax2.plot(Rx_coords,Rx_force,c='green',linewidth=2,label='F')
        
        ax2.plot(Rx_coords,Rx_force_constant,c='red',linewidth=2,label='κ')
        
        ax2.set_ylabel('F,κ',fontsize=15)
        
        ax2.set_ylim(min(Rx_force_constant)-5,max(Rx_force)+5)
        
        plt.legend(loc='upper right')
        
        plt.axvline(alpha,c='black',linestyle='--')
        
        plt.axvline(gamma,c='black',linestyle='--')
               
        plt.axhline(0,c='black',alpha=0.40)
        
        plt.axvline(0,c='pink',alpha=0.45)
        
        
        plt.savefig(filename,dpi=500)
        
        #Creating symmetrical plot boundaries
        if abs(Rx_coords[-1])<abs(Rx_coords[0]):
            plt.xlim(-Rx_coords[-1],Rx_coords[-1])
        else:
            plt.xlim(Rx_coords[0],-Rx_coords[0])
            
        
        
        plt.show()
        
        
        #print('Activation energy (kcal/mol) = ',Energy[BETA]-Energy[0])
        #print('Reaction energy (kcal/mol)= ',Energy[-1]-Energy[0])
        #print('Eact,2(kcal/mol)= ',Energy[BETA]-Energy[alpha_index])
        #print('Eact,1(kcal/mol)= ',Energy[alpha_index]-Energy[0])        
        
        return Energy[BETA]-Energy[0],Energy[-1]-Energy[0],TS_extent
    
    
class RFA_2files():
    
    def __init__(self,file1,file2):
        
        self.file1=file1
        self.file2=file2
        
        
    def get_number_of_lines(self):
        
        file1=self.file1
        
        file2=self.file2
        
        if os.path.isfile(file1):
            
            with open(file1,'r') as my_file:
            
                count1=len(my_file.readlines())
                
        else:
            
            print(file1,' Not Found ')
            
            
        
        if os.path.isfile(file2):
            
            with open(file2,'r') as my_file:
            
                count2=len(my_file.readlines())
                
        else:
            
            print(file2,' Not Found ')
            
        return count1,count2
    
    
    def get_Rx_energies_VS_coordinates_reverse_dir(self):
        
        file=self.file1
        
        nl=self.get_number_of_lines()[0]
        
        ref_word='Energies reported relative to the TS energy'
        
        pos=None
        
        list_rx_coords_ener=[]
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                    
                        pos=my_file.tell()-len(lc)
                        
                        break
                    
                my_file.seek(pos,0)
                
                for i in range(5):
                    my_file.readline()
                
                
                count=0
                while True:
                    
                    lc=my_file.readline()
                    
                    if '-----------------' not in lc:
                        
                        list_rx_coords_ener.append([])
                        
                        list_rx_coords_ener[count].append(eval(lc.split()[2]))
                        
                        list_rx_coords_ener[count].append(eval(lc.split()[1]))
                        
                        count+=1
                        
                        #print(lc)
                        
                    else:
                        
                        break
                        
                        
            
        else:
            
            print('File Not Found')
            
            
        return list_rx_coords_ener
    
    
    def get_Rx_energies_VS_coordinates_forw_dir(self):
        
        file=self.file2
        
        nl=self.get_number_of_lines()[1]
        
        ref_word='Energies reported relative to the TS energy'
        
        pos=None
        
        list_rx_coords_ener=[]
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                    
                        pos=my_file.tell()-len(lc)
                        
                        break
                    
                my_file.seek(pos,0)
                
                for i in range(5):
                    my_file.readline()
                
                
                count=0
                while True:
                    
                    lc=my_file.readline()
                    
                    if '-----------------' not in lc:
                        
                        list_rx_coords_ener.append([])
                        
                        list_rx_coords_ener[count].append(eval(lc.split()[2]))
                        
                        list_rx_coords_ener[count].append(eval(lc.split()[1]))
                        
                        count+=1
                        
                        #print(lc)
                        
                    else:
                        
                        break
                        
                        
            
        else:
            
            print('File Not Found')
            
            
        return list_rx_coords_ener
    
    def get_geometries_file1(self):
        
        file=self.file1
        
        nl=self.get_number_of_lines()[0]
        
        ref_word1='Input orientation:'
        
        ref_word2='CHANGE IN THE REACTION COORDINATE = '
        
        pos1=[]
        
        pos2=[]
        
        count=0
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word1 in lc:
                        
                        pos1.append(my_file.tell()-len(lc))
                        
                    
                    if ref_word2 in lc:
                        
                        pos2.append(my_file.tell()-len(lc))
                        
                        
                geom_pos=[]
                
                for r in pos2:
                    
                    for s in pos1:
                        
                        if s<r:
                            
                            match=s
                        else:
                            
                            break
                                              
                    geom_pos.append(match)                           
                
                        
                geom_pos.insert(0,pos1[0])
                                
                irc_geoms=[]
                
                count_geoms=0
                
                for pos in geom_pos:
                    
                    irc_geoms.append([])
                    
                    my_file.seek(pos,0)
                    
                    for i in range(5):
                        
                        my_file.readline()
                        
                    count_atoms=0
                    
                    while True:
                        
                        lc=my_file.readline()
                    
                        if '------' not in lc:
                            
                            irc_geoms[count_geoms].append([])
                            
                            irc_geoms[count_geoms][count_atoms].append(dict_atoms[int(lc.split()[1])])
                            irc_geoms[count_geoms][count_atoms].append(eval(lc.split()[3]))
                            irc_geoms[count_geoms][count_atoms].append(eval(lc.split()[4]))
                            irc_geoms[count_geoms][count_atoms].append(eval(lc.split()[5]))
                            
                            count_atoms+=1
                            
                        else:
                            
                            break
                    
                    
                    count_geoms+=1
            
        else:
            
            print(file,'not found!')
            
        return irc_geoms
    
    
    def get_geometries_file2(self):
        
        file=self.file2
        
        nl=self.get_number_of_lines()[1]
        
        ref_word1='Input orientation:'
        
        ref_word2='CHANGE IN THE REACTION COORDINATE = '
        
        pos1=[]
        
        pos2=[]
        
        count=0
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word1 in lc:
                        
                        pos1.append(my_file.tell()-len(lc))
                        
                    
                    if ref_word2 in lc:
                        
                        pos2.append(my_file.tell()-len(lc))
                        
                        
                geom_pos=[]
                
                for r in pos2:
                    
                    for s in pos1:
                        
                        if s<r:
                            
                            match=s
                        else:
                            
                            break
                                              
                    geom_pos.append(match)                           
                
                        
                geom_pos.insert(0,pos1[0])
                          
                
                irc_geoms=[]
                
                count_geoms=0
                
                for pos in geom_pos:
                    
                    irc_geoms.append([])
                    
                    my_file.seek(pos,0)
                    
                    for i in range(5):
                        
                        my_file.readline()
                        
                    
                       
                    count_atoms=0
                    
                    while True:
                        
                        lc=my_file.readline()
                    
                        if '------' not in lc:
                            
                            irc_geoms[count_geoms].append([])
                            
                            irc_geoms[count_geoms][count_atoms].append(dict_atoms[int(lc.split()[1])])
                            irc_geoms[count_geoms][count_atoms].append(eval(lc.split()[3]))
                            irc_geoms[count_geoms][count_atoms].append(eval(lc.split()[4]))
                            irc_geoms[count_geoms][count_atoms].append(eval(lc.split()[5]))
                            
                            count_atoms+=1
                            
                        else:
                            
                            break
                    
                    
                    count_geoms+=1
                
            
        else:
            
            print(file,'not found!')
            
            
        
        return irc_geoms
        
    
    def construct_geoms_string(self):
        
        geom1=self.get_geometries_file1()
        geom2=self.get_geometries_file2()
        
        dir1=self.get_Rx_energies_VS_coordinates_forw_dir()
        dir2=self.get_Rx_energies_VS_coordinates_reverse_dir()
        
                
        gap1=abs(dir1[-1][1]-dir1[0][1])
        gap2=abs(dir2[-1][1]-dir2[0][1])
        
        
        
        if gap2>gap1:
            
            irc_path=None
            
            
            if dir1[0][0]>=0 and dir1[1][0]>0:
                geom1_new=self.reverse_list(geom1)
                
            else:
                geom1_new=geom1
                
            if dir2[0][0]<0 and dir2[1][0]<0:                
                geom2_new=self.reverse_list(geom2)
            else:
                geom2_new=geom2
                
            rev=geom1_new
            forw=geom2_new
                
        else:
            
                        
            if dir2[0][0]>=0 and dir2[1][0]>0:
                geom2_new=self.reverse_list(geom2)
            else:
                geom2_new=geom2
                
            if dir1[0][0]<0 and dir1[1][0]<0:
                geom1_new=self.reverse_list(geom1)
            else:
                geom1_new=geom1
                
            rev=geom2_new
            forw=geom1_new
            
        irc_path=rev+forw[1:]
                
        return irc_path            
                
    
    def reverse_list(self,my_list):
        
        new_list=[]
        
        for i in range(len(my_list)):
            
            new_list.append(my_list[len(my_list)-1-i])
            
        return new_list
        
    def get_Rx_energies_VS_coordinates(self):
        
        dir1=self.get_Rx_energies_VS_coordinates_forw_dir()
        dir2=self.get_Rx_energies_VS_coordinates_reverse_dir()
        
                
        gap1=abs(dir1[-1][1]-dir1[0][1])
        gap2=abs(dir2[-1][1]-dir2[0][1])
        
        
        
        if gap2>gap1:
            
            dir1_new=[]
            dir2_new=[]
            
            if dir1[0][0]>=0 and dir1[1][0]>0:
                
                count=0                
                for k in range(len(dir1)):
                    dir1_new.append([])
                    
                    if dir1[k][0]!=0:
                    
                        dir1_new[count].append(-dir1[k][0])
                        dir1_new[count].append(dir1[k][1])
                    else:
                        dir1_new[count].append(dir1[k][0])
                        dir1_new[count].append(dir1[k][1])
                        
                    count+=1
                    
                    
                dir1_new=self.reverse_list(dir1_new)
                #dir1_new=dir1_new.reverse()
                    
            else:
                
                dir1_new=dir1             
            
            
            if dir2[0][0]<0 and dir2[1][0]<0:

                count=0
                for r in range(len(dir2)):
                    dir2_new.append([])
                    
                    if dir2[r][0]!=0:
                    
                        dir2_new[count].append(-dir2[r][0])
                        dir2_new[count].append(dir2[r][1])
                    else:
                        dir2_new[count].append(dir2[r][0])
                        dir2_new[count].append(dir2[r][1])
                    
                    
                    count+=1
                    
                dir2_new=self.reverse_list(dir2_new)
                #dir2_new=dir2_new.reverse()
                
            else:
                
                dir2_new=dir2
                
                
            
            rev=dir1_new
            forw=dir2_new
            
            
        else:
            dir1_new=[]
            dir2_new=[]
            
            if dir2[0][0]>=0 and dir2[1][0]>0:

                count=0
                for r in range(len(dir2)):
                    dir2_new.append([])
                    
                    if dir2[r][0]!=0:
                    
                        dir2_new[count].append(-dir2[r][0])
                        dir2_new[count].append(dir2[r][1])
                    else:
                        dir2_new[count].append(dir2[r][0])
                        dir2_new[count].append(dir2[r][1])
                    
                    count+=1
                    
                dir2_new=self.reverse_list(dir2_new)
                
            else:
                
                dir2_new=dir2                    
              
            
            if dir1[0][0]<0 and dir1[1][0]<0:
                
                count=0                
                for k in range(len(dir1)):
                    dir1_new.append([])
                    
                    if dir1[k][0]!=0:
                    
                        dir1_new[count].append(-dir1[k][0])
                        dir1_new[count].append(dir1[k][1])
                    else:
                        dir1_new[count].append(dir1[k][0])
                        dir1_new[count].append(dir1[k][1])
                        
                    count+=1
                    
                dir1_new=self.reverse_list(dir1_new)
                
            else:
                
                dir1_new=dir1
            
            rev=dir2_new
            forw=dir1_new
                
        irc=rev+forw[1:]
        
        return irc
    
    def get_Rx_coords(self):
        
        E_Rc=self.get_Rx_energies_VS_coordinates()
        
        Rc=[]
        
        for i in range(len(E_Rc)):
            
            Rc.append(E_Rc[i][0])
            
        return Rc
    
    def get_Energy(self):
        
        E_Rx_coords=self.get_Rx_energies_VS_coordinates()
        
        E=[]
        
        for i in range(len(E_Rx_coords)):
            
            E.append(E_Rx_coords[i][1])
            
        E=[(x-E[0])*627.509 for x in E]   
        
        return E 
    
    def get_Rx_force(self):
        
        E=self.get_Energy()
        
        Rc=self.get_Rx_coords()
        
        centr_der=[]
        
        first_term=(E[1]-E[0])/(Rc[1]-Rc[0])
        
        centr_der.append(first_term)
        
        for i in range(1,len(E)-1):
            
            forw_der=(E[i+1]-E[i])/(Rc[i+1]-Rc[i])
            
            backw_der=(E[i]-E[i-1])/(Rc[i]-Rc[i-1])
            
            centr_der.append((forw_der+backw_der)/2)
            
        centr_der.append((E[-2]-E[-1])/(Rc[-2]-Rc[-1]))
        
        force=[-x for x in centr_der]
        
        return force
    
    def get_min_max_Rx_force(self):
        
        force=self.get_Rx_force()
        
        Rc=self.get_Rx_coords()
        
        global_min=Rc[0]
        
        global_max=Rc[0]
        
        test=force[0]
        
        for i in range(1,len(force)):
            
            if force[i]<test:
                
                global_min=Rc[i]
                
                test=force[i]
                
        for j in range(1,len(force)):
                
            if force[j]>test:
                
                global_max=Rc[j]
                
                test=force[j]
                
        return global_min,global_max
    
    
    def get_Rx_force_constant(self):
        
        force=self.get_Rx_force()
        
        Rc=self.get_Rx_coords()
        
        centr_der=[(force[1]-force[0])/(Rc[1]-Rc[0])]
        
        for i in range(1,len(force)-1):
            
            forw_der=(force[i+1]-force[i])/(Rc[i+1]-Rc[i])
            
            backw_der=(force[i]-force[i-1])/(Rc[i]-Rc[i-1])
            
            centr_der.append((forw_der+backw_der)/2)
            
        centr_der.append((force[-2]-force[-1])/(Rc[-2]-Rc[-1]))
            
        force_constant=[-x for x in centr_der]
        
        return force_constant
    
        
    def get_TS_region_extend(self):
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]
        
        return round(abs(gamma_limit-alpha_limit),2)
    
    def dec_activ_energy(self):

        Alpha_index=(self.get_Rx_coords()).index(self.get_min_max_Rx_force()[0])
        Gamma_index=(self.get_Rx_coords()).index(self.get_min_max_Rx_force()[1])
        Energy=self.get_Energy()[0]
        E_TS=max(Energy)
        E_Reactant=Energy[0]
        E_Alpha=Energy[Alpha_index]
        E_Gamma=Energy[Gamma_index]
        E_Prod=Energy[-1]
        
        return round(E_Alpha-E_Reactant,2),round(E_TS-E_Alpha,2),round(E_Gamma -E_TS,2),round(E_Prod-E_Gamma,2)  
    
    
    def plot(self,output_dir):
        
        import matplotlib.pyplot as plt        
        import pandas as pd
        import numpy as np
        
        file_tail=os.path.split(self.file1)[1]
                
        filename=os.path.join(output_dir,file_tail[3:-4])
        
        title=(os.path.split(self.file1)[1])[5:-6]
        
        alpha=self.get_min_max_Rx_force()[0]
        
        gamma=self.get_min_max_Rx_force()[1]
        
        alpha_index=(self.get_Rx_coords()).index(self.get_min_max_Rx_force()[0])
        
        BETA=(self.get_Energy()).index(max(self.get_Energy()))
        
        GAMMA=(self.get_Rx_coords()).index(self.get_min_max_Rx_force()[1])
        
        TS_extent=self.get_TS_region_extend()
        
        Energy=self.get_Energy()
        
        Rx_force=self.get_Rx_force()
        
        Rx_force_constant=self.get_Rx_force_constant()
        
        Rx_coords=self.get_Rx_coords()        
        
        plt.figure(figsize=(6,3),dpi=3000)
        
        plt.title(file_tail[3:])
        
        fig,ax=plt.subplots()        
        
        plt.ylabel('Energy (kcal/mol)',fontsize=13)
        
        ax.plot(Rx_coords,Energy,c='blue',linewidth=2,label='E')
                
        ax.set_xlabel('Reaction coordinate (ξ)',labelpad=8)
        
        ax.set_ylim(min(Energy)-5,max(Energy)+5)
        
        plt.legend(loc='upper left')
        
        plt.axhline(0,c='black',alpha=0.40)
        
        plt.axvspan(min(Rx_coords),alpha,alpha=0.3,color='silver')
        
        plt.axvspan(gamma,max(Rx_coords),alpha=0.3,color='silver')
        
        ax.text(alpha-1.5,min(Energy)-3.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        ax.text(0.2,min(Energy)-3.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        ax.text(gamma+0.5,min(Energy)-3.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})

        ax2=ax.twinx()
        
        ax2.plot(Rx_coords,Rx_force,c='green',linewidth=2,label='F')
        
        ax2.plot(Rx_coords,Rx_force_constant,c='red',linewidth=2,label='κ')
        
        ax2.set_ylabel('F,κ',fontsize=13)
        
        ax2.set_ylim(min(Rx_force_constant)-8,max(Rx_force)+5)
        
        plt.legend(loc='upper right')
        
        plt.axvline(alpha,c='black',linestyle='--')
        
        plt.axvline(gamma,c='black',linestyle='--')
        
        plt.axhline(0,c='black',alpha=0.40)
        
        plt.axvline(0,c='pink',alpha=0.60)        
        
        plt.figure(figsize=(6,3),dpi=3000)     
                
        plt.savefig(filename,dpi=1500)
        
        plt.show()        
        
        #print('Activation energy (kcal/mol) = ',Energy[BETA]-Energy[0])
        #print('Reaction energy (kcal/mol)= ',Energy[-1]-Energy[0])
        #print('Eact,2(kcal/mol)= ',Energy[BETA]-Energy[alpha_index])
        #print('Eact,1(kcal/mol)= ',Energy[alpha_index]-Energy[0])
        #print('TS_extent (ξ)= ',TS_extent)
        
        return Energy[BETA]-Energy[0],Energy[-1]-Energy[0],Energy[BETA]-Energy[alpha_index],Energy[alpha_index]-Energy[0],TS_extent
    
class MultiPotentialEnergy2F():
    
    def __init__(self,list_tuples):
        
        self.list_tuples=list_tuples
        
        
    def get_PE_and_IRC_limits(self):
        
        files=self.list_tuples
        
        PEs=[]
        
        Rx_coords=[]
        
        IRC_sections=[]
        
        i=0
        
        for reaction in files:
            
            PEs.append([])
            
            Rx_coords.append([])
            
            IRC_sections.append([])
            
            current_reaction_files=[]
            
            for file in reaction:
                
                current_reaction_files.append(file)
                
            PEs[i].append(RFA_2files(current_reaction_files[0],current_reaction_files[1]).get_Energy())
            
            Rx_coords[i].append(RFA_2files(current_reaction_files[0],current_reaction_files[1]).get_Rx_coords())        
            
            IRC_sections[i].append(RFA_2files(current_reaction_files[0],current_reaction_files[1]).get_min_max_Rx_force())   
                
                
            i+=1
            
        return PEs,Rx_coords
    
    def plot(self):
        
        import matplotlib.pyplot as plt
        
        colors=['b','r','y','g','c','m','pink','silver']
        
        PEs=self.get_PE_and_IRC_limits()[0]
        
        Rc=self.get_PE_and_IRC_limits()[1]
        
        IRC_sections=self.get_PE_and_IRC_limits()[2]
        
        plt.figure(figsize=(60,20))
        
        plt.ylabel('Potential Energy (kcal/mol)',fontsize=8)
        
        plt.xlabel('IRC (ξ)',labelpad=8)
        
        plt.axhline(0,c='black',alpha=0.40)
        
        for k in range(len(PEs)):
            
            Energy=PEs[k]
            
            Rx_coords=Rc[k]
            
            alpha=IRC_sections[k][0]
            
            gamma=IRC_sections[k][1]
            
            plt.plot(Rx_coords,Energy,linewidth=1,label='Reaction'+str(k+1))
            
            plt.axvline(alpha,c=colors[k],linestyle='--',alpha=0.4)
        
            plt.axvline(gamma,c=colors[k],linestyle='--',alpha=0.4)
            
                      
        plt.legend(loc='best')
        
        
class MultiReactionForce2F():
    
    def __init__(self,list_tuples):
        
        self.list_tuples=list_tuples
        
        
    def get_PE_and_IRC_limits(self):
        
        files=self.list_tuples
        
        RFs=[]
        
        Rx_coords=[]
        
        IRC_sections=[]
        
        i=0
        
        for reaction in files:
            
            RFs.append([])
            
            Rx_coords.append([])
            
            IRC_sections.append([])
            
            current_reaction_files=[]
            
            for file in reaction:
                
                current_reaction_files.append(file)
                
            RFs[i].append(RFA_2files(current_reaction_files[0],current_reaction_files[1]).get_Rx_force())
            
            Rx_coords[i].append(RFA_2files(current_reaction_files[0],current_reaction_files[1]).get_Rx_coords())        
            
            IRC_sections[i].append(RFA_2files(current_reaction_files[0],current_reaction_files[1]).get_min_max_Rx_force())   
                
                
            i+=1
            
        return RFs,Rx_coords
    
    def plot(self):
        
        import matplotlib.pyplot as plt
        
        colors=['b','r','y','g','c','m','pink','silver']
        
        RFs=self.get_PE_and_IRC_limits()[0]
        
        Rc=self.get_PE_and_IRC_limits()[1]
        
        IRC_sections=self.get_PE_and_IRC_limits()[2]
        
        plt.figure(figsize=(60,20))
        
        plt.ylabel('Reaction Force (kcal/mol.ξ)',fontsize=8)
        
        plt.xlabel('IRC (ξ)',labelpad=8)
        
        plt.axhline(0,c='black',alpha=0.40)
        
        for k in range(len(RFs)):
            
            RF=RFs[k]
            
            Rx_coords=Rc[k]
            
            alpha=IRC_sections[k][0]
            
            gamma=IRC_sections[k][1]
            
            plt.plot(Rx_coords,RF,linewidth=1,label='Reaction'+str(k+1))
            
            plt.axvline(alpha,c=colors[k],linestyle='--',alpha=0.4)
        
            plt.axvline(gamma,c=colors[k],linestyle='--',alpha=0.4)
            
                      
        plt.legend(loc='best')
        
        
class MultiReactionForceConstant2F():
    
    def __init__(self,list_tuples):
        
        self.list_tuples=list_tuples
        
        
    def get_PE_and_IRC_limits(self):
        
        files=self.list_tuples
        
        RFCs=[]
        
        Rx_coords=[]
        
        IRC_sections=[]
        
        i=0
        
        for reaction in files:
            
            RFCs.append([])
            
            Rx_coords.append([])
            
            IRC_sections.append([])
            
            current_reaction_files=[]
            
            for file in reaction:
                
                current_reaction_files.append(file)
                
            RFCs[i].append(RFA_2files(current_reaction_files[0],current_reaction_files[1]).get_Rx_force_constant())
            
            Rx_coords[i].append(RFA_2files(current_reaction_files[0],current_reaction_files[1]).get_Rx_coords())        
            
            IRC_sections[i].append(RFA_2files(current_reaction_files[0],current_reaction_files[1]).get_min_max_Rx_force())   
                
                
            i+=1
            
        return RFCs,Rx_coords
    
    def plot(self):
        
        import matplotlib.pyplot as plt
        
        colors=['b','r','y','g','c','m','pink','silver']
        
        RFCs=self.get_PE_and_IRC_limits()[0]
        
        Rc=self.get_PE_and_IRC_limits()[1]
        
        IRC_sections=self.get_PE_and_IRC_limits()[2]
        
        plt.figure(figsize=(60,20))
        
        plt.ylabel('Reaction Force constant (kcal/mol.ξ^2)',fontsize=8)
        
        plt.xlabel('IRC (ξ)',labelpad=8)
        
        plt.axhline(0,c='black',alpha=0.40)
        
        for k in range(len(RFCs)):
            
            RFC=RFCs[k]
            
            Rx_coords=Rc[k]
            
            alpha=IRC_sections[k][0]
            
            gamma=IRC_sections[k][1]
            
            plt.plot(Rx_coords,RFC,linewidth=1,label='Reaction'+str(k+1))
            
            plt.axvline(alpha,c=colors[k],linestyle='--',alpha=0.4)
        
            plt.axvline(gamma,c=colors[k],linestyle='--',alpha=0.4)
            
                      
        plt.legend(loc='best')
        
        
class MultiPotentialEnergy1F():
    
    def __init__(self,list_files):
        
        self.list_files=list_files
        
        
    def get_PE_and_IRC_limits(self):
        
        files=self.list_files
        
        PEs=[]
        
        Rx_coords=[]
        
        IRC_sections=[]
        
        i=0
        
        for file in files:
            
            instance=RFA_1file(file)
                               
            PEs.append(instance.get_Energy()[0])
            
            Rx_coords.append(instance.get_Rx_coords())        
            
            IRC_sections.append(instance.get_min_max_Rx_force())   
               
            
        return PEs,Rx_coords,IRC_sections
    
    def plot(self,name,figs):
        
        import matplotlib.pyplot as plt
        
        colors=['b','r','y','g','c','m','pink','silver','#141da2',
               '#9c21d1','#d12182']
        
        Rx_IDs=[os.path.split(file)[1][3:-4] for file in self.list_files]
        
        PEs=self.get_PE_and_IRC_limits()[0]
        
        Rc=self.get_PE_and_IRC_limits()[1]
        
        IRC_sections=self.get_PE_and_IRC_limits()[2]
        
        plt.figure(figsize=(figs[0],figs[1]))
        
        plt.ylabel('Potential Energy (kcal/mol)',fontsize=8)
        
        plt.xlabel('Reaction coordinate (ξ)',labelpad=8)
        
        plt.axhline(0,c='black',alpha=0.40)
        
        plt.axvline(0,c='purple',alpha=0.4)
        
        for k in range(len(PEs)):
            
            Energy=PEs[k]
            
            Rx_coords=Rc[k]
            
            alpha=IRC_sections[k][0]
            
            gamma=IRC_sections[k][1]
            
            plt.plot(Rx_coords,Energy,c=colors[k],linewidth=1,label=Rx_IDs[k])
            
            plt.axvline(alpha,c=colors[k],linestyle='--')
        
            plt.axvline(gamma,c=colors[k],linestyle='--')
            
                      
        plt.legend(loc='best')
        plt.savefig(name,dpi=2000)
        
class MultiReactionForce1F():
    
    def __init__(self,list_files):
        
        self.list_files=list_files
        
        
    def get_PE_and_IRC_limits(self):
        
        files=self.list_files
        
        RFs=[]
        
        Rx_coords=[]
        
        IRC_sections=[]
        
        for file in files:
            
            instance=RFA_1file(file)
                               
            RFs.append(instance.get_Rx_force()[0])
            
            Rx_coords.append(instance.get_Rx_coords())        
            
            IRC_sections.append(instance.get_min_max_Rx_force())   
                
                    
        return RFs,Rx_coords,IRC_sections
    
    def plot(self,name,figs):
        
        import matplotlib.pyplot as plt
        
        colors=['b','r','y','g','c','m','pink','silver','#141da2',
               '#9c21d1','#d12182']
        
        Rx_IDs=[os.path.split(file)[1][3:-4] for file in self.list_files]
        
        RFs=self.get_PE_and_IRC_limits()[0]
        
        Rc=self.get_PE_and_IRC_limits()[1]
        
        IRC_sections=self.get_PE_and_IRC_limits()[2]
        
        plt.figure(figsize=(figs[0],figs[1]))
        
        plt.ylabel('Reaction Force (kcal/mol.ξ)',fontsize=8)
        
        plt.xlabel('Reaction coordinate (ξ)',labelpad=8)
        
        plt.axhline(0,c='black',alpha=0.40)
        
        plt.axvline(0,c='purple',alpha=0.4)
       
        for k in range(len(RFs)):
            
            RF=RFs[k]
            
            Rx_coords=Rc[k]
            
            alpha=IRC_sections[k][0]
            
            gamma=IRC_sections[k][1]
            
            plt.plot(Rx_coords,RF,c=colors[k],linewidth=1,label=Rx_IDs[k])
            
            plt.axvline(alpha,c=colors[k],linestyle='--')
        
            plt.axvline(gamma,c=colors[k],linestyle='--')
            
                        
                      
        plt.legend(loc='best')
        plt.savefig(name,dpi=2000)
        
        
class MultiReactionForceConstant1F():
    
    def __init__(self,list_files):
        
        self.list_files=list_files
        
        
    def get_PE_and_IRC_limits(self):
        
        files=self.list_files
        
        RFCs=[]
        
        Rx_coords=[]
        
        IRC_sections=[]
        
        for file in files:
            
            instance=RFA_1file(file)
                               
            RFCs.append(instance.get_Rx_force_constant()[0])
            
            Rx_coords.append(instance.get_Rx_coords())        
            
            IRC_sections.append(instance.get_min_max_Rx_force()) 
            
            
        return RFCs,Rx_coords,IRC_sections
    
    def plot(self,name,figs):
        
        import matplotlib.pyplot as plt
        
        colors=['b','r','y','g','c','m','pink','silver','#141da2',
               '#9c21d1','#d12182']
        
        Rx_IDs=[os.path.split(file)[1][3:-4] for file in self.list_files]
        
        RFCs=self.get_PE_and_IRC_limits()[0]
        
        Rc=self.get_PE_and_IRC_limits()[1]
        
        IRC_sections=self.get_PE_and_IRC_limits()[2]
        
        plt.figure(figsize=(figs[0],figs[1]),dpi=2500)
        
        plt.ylabel('Reaction Force constant (kcal/mol.$ξ^2$)',fontsize=8)
        
        plt.xlabel('Reaction coordinate (ξ)',labelpad=8)
        
        plt.axhline(0,c='black',alpha=0.40)
        
        plt.xlim(-16,16)
        
        plt.axvline(0,c='purple',alpha=0.4)
                
        for k in range(len(RFCs)):
            
            RFC=RFCs[k]
            
            Rx_coords=Rc[k]
            
            alpha=IRC_sections[k][0]
            
            gamma=IRC_sections[k][1]
            
            plt.plot(Rx_coords,RFC,c=colors[k],linewidth=1,label=Rx_IDs[k])
            
            plt.axvline(alpha,c=colors[k],linestyle='--')
        
            plt.axvline(gamma,c=colors[k],linestyle='--')
            
                      
        plt.legend(loc='best',fontsize=8)
        plt.savefig(name,dpi=2000)
        
        
class create_IRC_input():
    
    def __init__(self,file):
        
        self.filename=file
    
    def get_number_of_lines(self):
        
        ts=open(self.filename)
        
        nl=len(ts.readlines())
        
        ts.close()
        
        return  nl
    
    def to_decimal(self,x,precision:int=10):
        
        return format(x,f".{precision}f").lstrip().rstrip('0')
    
    def check_convergence(self):
        
        file=self.filename
        
        ref_word1="Error termination request"
        
        ref_word2='Error termination via'
        
        ref_word3=' Normal termination'
        
        nl=self.get_number_of_lines()
        
        result='Still running'
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
            
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word1 in lc or ref_word2 in lc:
                        
                        result='Failed'
                        
                        break
                    
                    if ref_word3 in lc:
                        
                        result='Successful'
        else:
            
            print('File not Found')
            
         
        return result
                
                
    
    def locate_ts(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word1='SCF Done'
        
        ref_word2='Standard orientation'
        
        ref_word3='Optimization completed'
        
        ref_word4='!   Optimized Parameters   !'
        
        pos1_scf=[]
        
        pos2_stand_or=[]
        
        pos3_opt_compl=[]
        
        pos4_opt_param=[]
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word1 in lc:
                        
                        end_line_pos=my_file.tell()
                        
                        pos1_scf.append(end_line_pos-len(lc))
                        
                        
                    elif ref_word2 in lc:
                        
                        end_line_pos=my_file.tell()
                        
                        pos2_stand_or.append(end_line_pos-len(lc))
                        
                    elif ref_word3 in lc:
                        
                        end_line_pos=my_file.tell()
                        
                        pos3_opt_compl.append(end_line_pos-len(lc))
                        
                    elif ref_word4 in lc:
                        
                        end_line_pos=my_file.tell()
                        
                        pos4_opt_param.append(end_line_pos-len(lc))
                    else:
                        
                        None
            
            
        else:
            
            print('File Not Found')
            
        #Let us locate different configurations
        
        scf_done_line_pos=[]
        
        stand_orient_line_pos=[]
        
        opt_param_line_pos=[]
        
        for pos3 in pos3_opt_compl:
            
            scf_macth=0
            
            stand_or_match=0
            
            for pos1 in pos1_scf:
                             
                if pos1<pos3:
                    
                    scf_match=pos1
                    
            scf_done_line_pos=scf_match-2*len('SCF Done')
            
            for pos2 in pos2_stand_or:

                
                if pos2<pos3:
                    
                    stand_or_match=pos2
                    
                
            stand_orient_line_pos=stand_or_match-len('Standard orientation')
            
            for pos4 in pos4_opt_param:
                
                if pos4>pos3 and (pos4-pos3)<250:
                    
                    opt_param_match=pos4
                
            opt_param_line_pos=opt_param_match-len('!   Optimized Parameters   !')
            
                    
            
        return scf_done_line_pos,stand_orient_line_pos,opt_param_line_pos
    
    def locate_last_reached_geom(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Standard orientation'
        
        ref_word2='Charge =  0 Multiplicity = 1'
        
        pos=[0]
        
        pos2=[0]
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                    
                        pos.append(my_file.tell()-len(lc))
                        
                    if ref_word2 in lc:
                        
                        pos2.append(my_file.tell()-len(lc))
                        
        else:
            
            print(f'{file} not found')
        
        
            
        return pos[-1],pos2[-1]

    def get_input_geom(self):
        
        file=self.filename
        nl=self.get_number_of_lines()
        ref_word='Charge =  0 Multiplicity = 1'
        geom=[]
        pos=None
        with open(file,'r') as my_file:
            
            for i in range(nl):
                lc=my_file.readline()
                if ref_word in lc:
                    pos=my_file.tell()-len(lc)
                    break
            my_file.seek(pos,0)
            my_file.readline()
            count=0
            while True:
                geom.append([])
                
                lc=my_file.readline()
                
                if len(lc.split())>5:
                    
                    break
                else:
                    for i in range(len(lc.split())):
                        
                        geom[count].append(lc.split()[i])
                count+=1
                
        return geom[:-1]
                    
    def get_ts_geometry(self):
        
        #This function locate geometries of different stationary points 
        
        file=self.filename
        
        if self.check_convergence()=='Successful':
        
            stand_or_pos=self.locate_ts()[1]
            
        else:
            
            stand_or_pos=self.locate_last_reached_geom()[0]
            
        
        ts_geom=[]
        
        
        with open(file,'r') as my_file:
                
            my_file.seek(stand_or_pos,0)
            
            for i in range(6):
                
                my_file.readline()
                
                            
            lc=my_file.readline()
                
            line_count=0
                
            while '----------' not in lc:
                
                ts_geom.append([])
                    
                    
                for i in range(len(lc.split())):
                    
                    value=eval(lc.split()[i])
                    
                    ts_geom[line_count].append(value)
                        
                line_count+=1
                    
                lc=my_file.readline()
                
                
                
                
        return ts_geom
    
    
    def Df_ts_geometry(self):
        
        import pandas as pd
        
        
        if self.check_convergence()=='Successful' or self.locate_last_reached_geom()[0]!=0:
        
            ts=self.get_ts_geometry()
            
            data=pd.DataFrame(ts,columns=['Center Number','Atomic Number','Atomic Type','X','Y','Z'])
            
        else:
            ts=self.get_input_geom()
            
            data=pd.DataFrame(ts,columns=['Atomic Symbol','X','Y','Z'])      
            
        
        return data
        
        
    def get_double_input(self,output_dir):
        
        file=self.filename
        
        TS=self.Df_ts_geometry()        
        
        file_tail=os.path.split(file)[1]
        
        new_filename1=os.path.join(output_dir,'IRC1'+file_tail[3:-3]+'gjf')
        
        with open(new_filename1,'w') as myfile:
            
            myfile.write('%nprocshared=8\n%mem=4GB\n# irc=(maxpoints=100,reverse,recorrect=test,calcall,stepsize=8) b3lyp/6-31g(d)\n\n')
            myfile.write('IRC\n\n0 1\n')
            
            if self.check_convergence()=='Successful':
                
                for index in TS.index:
                    
                    myfile.writelines([str(dict_atoms[TS['Atomic Number'][index]]),'\t\t',str(self.to_decimal(TS['X'][index])),'\t\t',str(self.to_decimal(TS['Y'][index])),'\t\t',str(self.to_decimal(TS['Z'][index])),'\n'])
                    
                myfile.write('\n')
                
            else:
                
                for index in TS.index:
                    
                    myfile.writelines([str([TS['Atomic Symbol'][index]]),'\t\t',str(self.to_decimal(TS['X'][index])),'t\t',str(self.to_decimal(TS['Y'][index])),'\)t\t',str(self.to_decimal(TS['Z'][index])),'\n'])
                    
                myfile.write('\n')
                
        
        new_filename2=os.path.join(output_dir,'IRC2'+file_tail[3:-3]+'gjf')
        
        with open(new_filename2,'w') as myfile:
            
            myfile.write('%nprocshared=8\n%mem=4GB\n# irc=(maxpoints=100,forward,recorrect=test,calcall,stepsize=8) b3lyp/6-31g(d)\n\n')
            myfile.write('IRC\n\n0 1\n')
            
            if self.check_convergence()=='Successful':
                
                for index in TS.index:
                    
                    myfile.writelines([str(dict_atoms[TS['Atomic Number'][index]]),'\t\t',str(self.to_decimal(TS['X'][index])),'\t\t',str(self.to_decimal(TS['Y'][index])),'\t\t',str(self.to_decimal(TS['Z'][index])),'\n'])
                    
                myfile.write('\n')
                
            else:
                
                for index in TS.index:
                    
                    myfile.writelines([str([TS['Atomic Symbol'][index]]),'\t\t',str(self.to_decimal(TS['X'][index])),'t\t',str(self.to_decimal(TS['Y'][index])),'\)t\t',str(self.to_decimal(TS['Z'][index])),'\n'])
                    
                myfile.write('\n')
                
                
    def get_unique_input(self,output_dir):
        
        file=self.filename
            
        TS=self.Df_ts_geometry()        
            
        file_tail=os.path.split(file)[1]
            
        new_filename=os.path.join(output_dir,'IRC'+file_tail[3:-3]+'gjf')
            
        with open(new_filename,'w') as myfile:
                
            myfile.write('%nprocshared=8\n%mem=4GB\n# irc=(maxpoints=100,recorrect=test,calcall,stepsize=8) b3lyp/6-31g(d)\n\n')
                
            myfile.write('IRC\n\n0 1\n')
                
            if self.check_convergence()=='Successful':
                    
                for index in TS.index:
                        
                    myfile.writelines([str(dict_atoms[TS['Atomic Number'][index]]),'\t\t',str(self.to_decimal(TS['X'][index])),'\t\t',str(self.to_decimal(TS['Y'][index])),'\t\t',str(self.to_decimal(TS['Z'][index])),'\n'])
                        
                myfile.write('\n')
                    
            else:
                    
                for index in TS.index:
                        
                    myfile.writelines([str([TS['Atomic Symbol'][index]]),'\t\t',str(self.to_decimal(TS['X'][index])),'t\t',str(self.to_decimal(TS['Y'][index])),'\)t\t',str(self.to_decimal(TS['Z'][index])),'\n'])
                        
                myfile.write('\n')
    
    
    
    
    