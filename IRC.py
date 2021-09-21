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



import os
import shutil
from tqdm import tqdm 

dict_atoms={
    
    1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C',
    
    7: 'N', 8:'O', 9:'F', 10:'Ne', 11:'Na', 12: 'Mg',
    
    13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl', 18: 'Ar',
    
    19: 'K', 20: 'Ca', 21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',
    
    27:'Co',28:'Ni',29:'Cu',30:'Zn',31:'Ga',32:'Ge',33:'As',34:'Se',
    
    35:'Br',36:'Kr',37:'Rb',38:'Sr',39:'Y',40:'Zr',41:'Nb',42:'Mo',43:'Tc',
    
    44:'Ru',45:'Rh',46:'Pd',47:'Ag',48:'Cd',49:'In',50:'Sn',51:'Sb',52:'Te',
    
    53:'I',54:'Xe'}



class IRC_1F():
    
    def __init__(self,filename):
        
        self.filename=filename
        
    def get_file_size(self,Mb=True):
        
        file=self.filename
        
        if os.path.isfile(file):
            
            if Mb==True:
                
                size=(os.stat(file).st_size)/(1024*1024)
                
            else:
                
                size=(os.stat(file).st_size)/(1024)
        
        else:                
            
            print('File Not found')
            
            
        return round(size,2)
    
    
    def get_number_of_lines(self):
        
        file=self.filename
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
            
                count=len(my_file.readlines())
                
        else:
            
            print('File Not Found ')
            
        return count
    
    def reverse(self,my_list):
        
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
    
    
    def get_geometries_file(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word1='Input orientation:'
        
        ref_word2='CHANGE IN THE REACTION COORDINATE = '
        
        ref_word3='Path Number:   1'
        
        ref_word4='Path Number:   2'
        
        ref_word5='Calculation of REVERSE path complete.'
        
        pos1=[]
        
        pos2=[]
        
        path1_length=0
        
        path2_length=0
        
        pos_path2=[]
        
        count=0
        
        irc_geoms=[]
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word1 in lc:
                        
                        pos1.append(my_file.tell()-len(lc))
                        
                    
                    if ref_word2 in lc:
                        
                        pos2.append(my_file.tell()-len(lc))
                        
                        
                    if ref_word3 in lc:
                        
                        path1_length+=1
                        
                    if ref_word4 in lc:
                        
                        path2_length+=1     
                        
                        pos_path2.append(my_file.tell()-len(lc))
                        
                    if ref_word5 in lc:
                        
                        reverse=my_file.tell()-len(lc)
                        
                    
                
                        
                geom_pos=[]
                
                geom_dir1=[]
                
                geom_dir2=[]
                
                for r in pos2:
                    
                    for s in pos1:
                        
                        if s<r:
                            
                            match=s
                        else:
                            
                            break
                                              
                    geom_pos.append(match)                           
                
                        
                geom_pos.insert(0,pos1[0])
                                
                
                
                count_geom1=0
                
                count_geom2=0
                
                count=0 
                
                for pos in geom_pos:
                    
                    if count<path1_length:
                    
                        geom_dir1.append([])
                        
                        my_file.seek(pos,0)
                        
                        for i in range(5):
                            
                            my_file.readline()
                            
                        count_atoms=0
                        
                        while True:
                            
                            lc=my_file.readline()
                        
                            if '------' not in lc:
                                
                                geom_dir1[count_geom1].append([])
                                
                                geom_dir1[count_geom1][count_atoms].append(dict_atoms[int(lc.split()[1])])
                                geom_dir1[count_geom1][count_atoms].append(eval(lc.split()[3]))
                                geom_dir1[count_geom1][count_atoms].append(eval(lc.split()[4]))
                                geom_dir1[count_geom1][count_atoms].append(eval(lc.split()[5]))
                                
                                count_atoms+=1
                                
                            else:
                                
                                break
                        
                        
                        count_geom1+=1
                        
                    else:
                        
                        geom_dir2.append([])
                        
                        my_file.seek(pos,0)
                        
                        for i in range(5):
                            
                            my_file.readline()
                            
                        count_atoms=0
                        
                        while True:
                            
                            lc=my_file.readline()
                        
                            if '------' not in lc:
                                
                                geom_dir2[count_geom2].append([])
                                
                                geom_dir2[count_geom2][count_atoms].append(dict_atoms[int(lc.split()[1])])
                                geom_dir2[count_geom2][count_atoms].append(eval(lc.split()[3]))
                                geom_dir2[count_geom2][count_atoms].append(eval(lc.split()[4]))
                                geom_dir2[count_geom2][count_atoms].append(eval(lc.split()[5]))
                                
                                count_atoms+=1
                                
                            else:
                                
                                break
                        
                        
                        count_geom2+=1
                        
                    count+=1
                
                if reverse>pos_path2[-1]:
                    
                    new_geom2=self.reverse(geom_dir2)
                    new_geom1=geom_dir1
                    irc_geoms=new_geom2+new_geom1
                    
                else:
                    
                    new_geom1=self.reverse(geom_dir1)
                    new_geom2=geom_dir2
                    irc_geoms=new_geom1+new_geom2
                    
            
        else:
            
            print(file,'not found!')
            
        return irc_geoms
    
    
    def get_atomic_force_components(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word1='Forces (Hartrees/Bohr)'
        
        ref_word2='CHANGE IN THE REACTION COORDINATE = '
        
        ref_word3='Path Number:   1'
        
        ref_word4='Path Number:   2'
        
        ref_word5='Calculation of REVERSE path complete.'
        
        pos1=[]
        
        pos2=[]
        
        path1_length=0
        
        path2_length=0
        
        pos_path2=[]
        
        count=0
        
        irc_force_matrix=[]
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word1 in lc:
                        
                        pos1.append(my_file.tell()-len(lc))
                        
                    
                    if ref_word2 in lc:
                        
                        pos2.append(my_file.tell()-len(lc))
                        
                        
                    if ref_word3 in lc:
                        
                        path1_length+=1
                        
                    if ref_word4 in lc:
                        
                        path2_length+=1     
                        
                        pos_path2.append(my_file.tell()-len(lc))
                        
                    if ref_word5 in lc:
                        
                        reverse=my_file.tell()-len(lc)
                        
                    
                
                        
                force_matrix_pos=[]
                
                force_matrix_dir1=[]
                
                force_matrix_dir2=[]
                
                for r in pos2:
                    
                    for s in pos1:
                        
                        if s<r:
                            
                            match=s
                        else:
                            
                            break
                                              
                    force_matrix_pos.append(match)                           
                
                        
                force_matrix_pos.insert(0,pos1[0])
                                
                
                
                count_fm1=0
                
                count_fm2=0
                
                count=0 
                
                for pos in force_matrix_pos:
                    
                    if count<path1_length:
                    
                        force_matrix_dir1.append([])
                        
                        my_file.seek(pos,0)
                        
                        for i in range(3):
                            
                            my_file.readline()
                            
                        count_atoms=0
                        
                        while True:
                            
                            lc=my_file.readline()
                        
                            if '------' not in lc:
                                
                                force_matrix_dir1[count_fm1].append([])
                                
                                force_matrix_dir1[count_fm1][count_atoms].append(eval(lc.split()[0]))
                                force_matrix_dir1[count_fm1][count_atoms].append(dict_atoms[int(lc.split()[1])])
                                force_matrix_dir1[count_fm1][count_atoms].append(eval(lc.split()[2]))
                                force_matrix_dir1[count_fm1][count_atoms].append(eval(lc.split()[3]))
                                force_matrix_dir1[count_fm1][count_atoms].append(eval(lc.split()[4]))
                                
                                count_atoms+=1
                                
                            else:
                                
                                break
                        
                        
                        count_fm1+=1
                        
                    else:
                        
                        force_matrix_dir2.append([])
                        
                        my_file.seek(pos,0)
                        
                        for i in range(3):
                            
                            my_file.readline()
                            
                        count_atoms=0
                        
                        while True:
                            
                            lc=my_file.readline()
                        
                            if '------' not in lc:
                                
                                force_matrix_dir2[count_fm2].append([])
                                
                                force_matrix_dir2[count_fm2][count_atoms].append(eval(lc.split()[0]))
                                force_matrix_dir2[count_fm2][count_atoms].append(dict_atoms[int(lc.split()[1])])
                                force_matrix_dir2[count_fm2][count_atoms].append(eval(lc.split()[2]))
                                force_matrix_dir2[count_fm2][count_atoms].append(eval(lc.split()[3]))
                                force_matrix_dir2[count_fm2][count_atoms].append(eval(lc.split()[4]))
                                
                                count_atoms+=1
                                
                            else:
                                
                                break
                        
                        
                        count_fm2+=1
                        
                    count+=1
                
                if reverse>pos_path2[-1]:
                    
                    new_force_matrix_2=self.reverse(force_matrix_dir2)
                    new_force_matrix_1=force_matrix_dir1
                    irc_force_matrix=new_force_matrix_2+new_force_matrix_1
                    
                else:
                    
                    new_force_matrix_1=self.reverse(force_matrix_dir1)
                    new_force_matrix_2=force_matrix_dir2
                    irc_force_matrix=new_force_matrix_1+new_force_matrix_2
                    
            
        else:
            
            print(file,'not found!')
            
        return irc_force_matrix
            
    def construct_geoms_string(self):
        
        geoms=self.get_geometries_file()
        
        direction=self.get_Rx_coords()[1]
        
        
        if direction=='normal':              
                            
            irc_geoms=geoms
            
        else:
            
            irc_geoms=self.reverse_list(geoms)              
    
        return irc_geoms
    
    def construct_irc_string_force_matrix(self):
        
        fms=self.get_atomic_force_components()
        
        direction=self.get_Rx_coords()[1]
    
        
        if direction=='normal':              
                            
            irc_fms=fms
            
        else:
            
            irc_fms=self.reverse_list(fms)              
    
        return irc_fms
    
    def Df_geoms_string(self):
        
        import pandas as pd
        
        
        geoms=self.construct_geoms_string()
        
        result=[]
        
        for geom in geoms:
            
            data=pd.DataFrame(geom,columns=['Atomic Symbol','X','Y','Z'])
            
            result.append(data)
    
        return result
    
    def Df_force_matrix_string(self):
        
        import pandas as pd
        
        
        fms=self.construct_irc_string_force_matrix()
        
        result=[]
        
        for fm in fms:
            
            
            data=pd.DataFrame(fm,columns=['Atomic Index','Atomic Symbol','Fx','Fy','Fz'])
            
            result.append(data)
    
        return result
    
    def reverse_list(self,my_list):
        
        new_list=[]
        
        for i in range(len(my_list)):
            
            new_list.append(my_list[len(my_list)-1-i])
            
        return new_list   
    
    def get_Rx_coords(self):
        
        E_Rc=self.get_Rx_energies_VS_coordinates()
        
        Rc=[]
        
        E=[]
        
        direction='normal'
        
        for i  in range(len(E_Rc)):
            
            Rc.append(E_Rc[i][0])
            
            E.append(E_Rc[i][1])
            
        if abs(E[0])>abs(E[-1]) and Rc[0]<Rc[-1]:
            
            Rc=[-x for x in self.reverse_list(Rc)]
            
            direction='reverse'
            
        if abs(E[0])>abs(E[-1]) and Rc[0]>Rc[-1]:
            
            Rc=self.reverse_list(Rc)
            
            direction='reverse'
        
                    
        return Rc,direction
    
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
            
            E=self.reverse(E)
            
            order='Reverse'
            
        
        return E,order
    
    def get_min_max_Rx_force(self):
        
        force=self.get_Rx_force()[0]
        
        Rc=self.get_Rx_force()[1]
        
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
    
    def plot_E_VS_Rx_coords(self):
        
        import matplotlib.pyplot as plt
        
        Rc=self.get_Rx_coords()

        Energy=self.get_Energy()[0]
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]

        plt.figure(figsize=(12,8))
        
        plt.ylim(min(Energy)-6,max(Energy)+6)
        
        plt.axvline(alpha_limit,0,color='black')
        
        plt.axvline(gamma_limit,0,color='black')
        
        plt.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plt.axhline(0)
        
        plt.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plt.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plt.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plt.scatter(Rc,Energy,c='blue',linewidths=0.5)
        
        plt.ylabel('Potential Energy, E(ξ)',fontsize=20)
        
        plt.xlabel('Reaction coordinate ξ',fontsize=20)
        
        #plt.text((Rc[-1]+gamma_limit+1)/2-2.8,max(Energy)+3,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow(Rc[-1]-gamma_limit,max(Energy)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow(Rc[-1]-gamma_limit,max(Energy)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        #plt.text((Rc[0]+alpha_limit)/2-3,max(Energy)+3,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(Energy)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(Energy)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        #plt.text((alpha_limit+gamma_limit)/2 - 1.5,max(Energy)+3,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(Energy)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(Energy)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.text(alpha_limit-0.5,min(Energy)-3.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        plt.text(0.4,min(Energy)-3.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plt.text(gamma_limit+0.5,min(Energy)-3.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
        
        
        plt.show()
        
        
    

    def get_Rx_force(self):
        
        E=self.get_Energy()[0]
        
        Rc=self.get_Rx_coords()[0]
        
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
        
        Rc=self.get_Rx_force()[1]
        
        centr_der=[(force[1]-force[0])/(Rc[1]-Rc[0])]
        
        for i in range(1,len(force)-1):
            
            forw_der=(force[i+1]-force[i])/(Rc[i+1]-Rc[i])
            
            backw_der=(force[i]-force[i-1])/(Rc[i]-Rc[i-1])
            
            centr_der.append((forw_der+backw_der)/2)
            
        centr_der.append((force[-2]-force[-1])/(Rc[-2]-Rc[-1]))
            
        force_constant=[-x for x in centr_der]
        
        return force_constant,Rc
 
    
    def plot_Rx_force(self):
        
        import matplotlib.pyplot as plt
        
        Rc=self.get_Rx_coords()
            
        forces=self.get_Rx_force()[0]
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]
        
        plt.figure(figsize=(12,8))
        
        plt.scatter(Rc,forces,c='blue')
        
        plt.ylim(min(forces)-3,max(forces)+3)
        
        plt.axvline(alpha_limit,0,color='black')
        
        plt.axvline(gamma_limit,0,color='black')
        
        #plt.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plt.axhline(0)
        
        plt.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plt.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plt.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plt.ylabel('Reaction force,F(ξ)',fontsize=20)
        
        #plt.xlabel('Reaction coordinate (ξ)',fontsize=20)
        
        plt.text((Rc[-1]+gamma_limit+1)/2-2.8,max(forces)+2,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow(Rc[-1]-gamma_limit,max(forces)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow(Rc[-1]-gamma_limit,max(forces)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        #plt.text((Rc[0]+alpha_limit)/2-3,max(forces)+2,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(forces)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(forces)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        #plt.text((alpha_limit+gamma_limit)/2 - 1.5,max(forces)+2,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(forces)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(forces)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.text(alpha_limit-0.5,min(forces)-2.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        #plt.text(0.4,min(forces)-2.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plt.text(gamma_limit+0.5,min(forces)-2.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
        #plt.savefig('fig.png')
        
        plt.show()
        
    
    
    def plot_Rx_force_constant(self):
        
        import matplotlib.pyplot as plt
        
        Rc=self.get_Rx_coords()
        
        force_constant=self.get_Rx_force_constant()[0]
        
        plt.figure(figsize=(12,8),dpi=600)
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]
        
        plt.ylim(min(force_constant)-3,max(force_constant)+3)
        
        plt.axvline(alpha_limit,0,color='black')
        
        plt.axvline(gamma_limit,0,color='black')
        
        #plt.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plt.axhline(0)
        
        plt.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plt.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plt.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plt.scatter(Rc,force_constant,c='blue')      
        
        plt.ylabel('Reactant force constant,κ(ξ)',fontsize=20)
        
        plt.xlabel('Reaction Coordinates (ξ)',fontsize=20)
        
        #plt.text((Rc[-1]+gamma_limit+1)/2-2.8,max(force_constant)+2,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow(Rc[-1]-gamma_limit,max(force_constant)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow(Rc[-1]-gamma_limit,max(force_constant)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        #plt.text((Rc[0]+alpha_limit)/2-3,max(force_constant)+2,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(force_constant)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(force_constant)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        #plt.text((alpha_limit+gamma_limit)/2 - 1.5,max(force_constant)+2,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(force_constant)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(force_constant)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.text(alpha_limit-0.5,min(force_constant)-2.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        #plt.text(0.4,min(force_constant)-2.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plt.text(gamma_limit+0.5,min(force_constant)-2.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
        
        plt.show()
    
    def plot_all(self):
        
        import matplotlib.pyplot as plt
        
        Rc=self.get_Rx_coords()
        
        Energy=self.get_Energy()[0]
        
        forces=self.get_Rx_force()[0]
        
        force_constant=self.get_Rx_force_constant()[0]
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]
        
        plt.figure(figsize=(32,45),dpi=(200))
        
        plot1=plt.subplot2grid((3,2),(0,0),colspan=2)
        
        plot2=plt.subplot2grid((3,2),(1,0),colspan=2)
        
        plot3=plt.subplot2grid((3,2),(2,0),colspan=2)
                
        #Plot 1
        
        plot1.set_ylim(min(Energy)-6,max(Energy)+6)
        
        plot1.axvline(alpha_limit,0,color='black')
        
        plot1.axvline(gamma_limit,0,color='black')
        
        plot1.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plot1.axhline(0)
        
        plot1.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plot1.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plot1.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plot1.scatter(Rc,Energy,c='blue',linewidths=0.5)
        
        plot1.set_ylabel('Potential Energy, E(ξ)',fontsize=20)
        
        plot1.set_xlabel('Reaction coordinate ξ',fontsize=20)
        
        #plot1.text(gamma_limit+1,max(Energy)+3,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot1.arrow(Rc[-1]-gamma_limit,max(Energy)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot1.arrow(Rc[-1]-gamma_limit,max(Energy)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        #plot1.text(Rc[0]+1,max(Energy)+3,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot1.arrow((Rc[0]+alpha_limit)/2,max(Energy)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot1.arrow((Rc[0]+alpha_limit)/2,max(Energy)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        #plot1.text((alpha_limit+gamma_limit)/2 - 1.5,max(Energy)+3,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot1.arrow((alpha_limit+gamma_limit)/2,max(Energy)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot1.arrow((alpha_limit+gamma_limit)/2,max(Energy)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot1.text(alpha_limit-0.5,min(Energy)-3.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        plot1.text(0.4,min(Energy)-3.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plot1.text(gamma_limit+0.5,min(Energy)-3.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
        
        #Plot2
        
        plot2.set_ylim(min(forces)-3,max(forces)+3)
        
        plot2.axvline(alpha_limit,0,color='black')
        
        plot2.axvline(gamma_limit,0,color='black')
        
        plot2.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plot2.axhline(0)
        
        plot2.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plot2.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plot2.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plot2.scatter(Rc,forces,color='blue')
        
        plot2.set_ylabel('Reaction force,F(ξ)',fontsize=20)
        
        plot2.set_xlabel('Reaction coordinate (ξ)',fontsize=20)
        
        #plot2.text(gamma_limit+1,max(forces)+2,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot2.arrow(Rc[-1]-gamma_limit,max(forces)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot2.arrow(Rc[-1]-gamma_limit,max(forces)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        #plot2.text(Rc[0]+1,max(forces)+2,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot2.arrow((Rc[0]+alpha_limit)/2,max(forces)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot2.arrow((Rc[0]+alpha_limit)/2,max(forces)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        #plot2.text((alpha_limit+gamma_limit)/2 - 1.5,max(forces)+2,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot2.arrow((alpha_limit+gamma_limit)/2,max(forces)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot2.arrow((alpha_limit+gamma_limit)/2,max(forces)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot2.text(alpha_limit-0.5,min(forces)-2.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        plot2.text(0.4,min(forces)-2.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plot2.text(gamma_limit+0.5,min(forces)-2.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
        
        #Plot 3
        
        plot3.set_ylim(min(force_constant)-3,max(force_constant)+3)
        
        plot3.axvline(alpha_limit,0,color='black')
        
        plot3.axvline(gamma_limit,0,color='black')
        
        plot3.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plot3.axhline(0)
        
        plot3.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plot3.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plot3.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plot3.scatter(Rc,force_constant,c='blue')      
        
        plot3.set_ylabel('Reactant force constant,κ(ξ)',fontsize=20)
        
        plot3.set_xlabel('Reaction Coordinates (ξ)',fontsize=20)
        
        #plot3.text(gamma_limit+1,max(force_constant)+2,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot3.arrow(Rc[-1]-gamma_limit,max(force_constant)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot3.arrow(Rc[-1]-gamma_limit,max(force_constant)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        #plot3.text(Rc[0]+1,max(force_constant)+2,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot3.arrow((Rc[0]+alpha_limit)/2,max(force_constant)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot3.arrow((Rc[0]+alpha_limit)/2,max(force_constant)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        #plot3.text((alpha_limit+gamma_limit)/2 - 1.5,max(force_constant)+2,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot3.arrow((alpha_limit+gamma_limit)/2,max(force_constant)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot3.arrow((alpha_limit+gamma_limit)/2,max(force_constant)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot3.text(alpha_limit-0.5,min(force_constant)-2.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        plot3.text(0.4,min(force_constant)-2.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plot3.text(gamma_limit+0.5,min(force_constant)-2.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
                
     
    def get_TS_region_extend(self):
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]
        
        return round(abs(gamma_limit-alpha_limit),2)
    
    def generate_irc_geom_spcalc_com_files(self,new_dir,level_of_theory):
        
        irc_path=self.construct_geoms_string()
        
        file_tail=os.path.split(self.filename)[1]
        
        count=1
        
        for geom in irc_path:
            
            filename=os.path.join(new_dir,file_tail[:-4]+'pt'+str(count)+'.gjf')
            
            with open(filename,'w') as my_file:
                
                my_file.write('%nprocshared=4\n')
                my_file.write('#p '+level_of_theory+' pop=nboread \n\n')
                my_file.write('SP\n\n0 1\n')
                
                for line in geom:
                    
                    my_file.writelines([str(line[0]),'\t\t',str(line[1]),'\t\t',str(line[2]),'\t\t',str(line[3]),'\n'])

                my_file.write('\n$nbo bndidx $end\n') 
                

            count+=1 

    def decompose_activation_energy(self):

        ALPHA=self.get_min_max_Rx_force()[0]
        BETA=0
        GAMMA=self.get_min_max_Rx_force()[1]
        Rc=self.get_Rx_coords()
        Energy=self.get_Energy()[0]
        Alpha_index=Rc.index(ALPHA)
        Beta_index=Rc.index(BETA)
        E_TS=Energy[Beta_index]
        E_Reactant=Energy[0]
        E_Alpha=Energy[Alpha_index]
        
        return round(E_TS-E_Reactant,2),round(E_Alpha-E_Reactant,2),round(E_TS-E_Alpha,2)
        
    
    
     
    def get_geoms_Dataframes(self):
        
        import pandas as pd
        
        irc_path=self.construct_geoms_string()
        IRC=[]
        for geom in irc_path:
            
            data=pd.DataFrame(geom,columns=['Atomic Symbol','X','Y','Z'])
            data['Atomic Label']=[x for x in range(1,len(data)+1)]
            IRC.append(data)
            
        return IRC 
    
    def get_reduced_geoms_Dataframes(self):  # These are H depleted geometries
    
        IRC=self.get_geoms_Dataframes()
        new_IRC=[]
        
        for geom in IRC:
            
            H_depleted_geom=geom[geom['Atomic Symbol']!='H']
            
            new_IRC.append(H_depleted_geom)
            
        return new_IRC
                
    
    def compute_distance(self,atom1,atom2):
        
        import math
        import numpy as np
        
        coords_zipper=zip(atom1,atom2)
        vect=[(r-s) for r,s in coords_zipper]        
        dist=np.linalg.norm(vect)
        
        return dist
    
    def extract_irc_react_and_prod_geom_com_files(self,new_dir,route_section='opt=(maxstep=1,notrustupdate,calcall) b3lyp/6-31G(d)'):
        
        IRC=self.get_geoms_Dataframes()
        Reactant=IRC[0]
        Product=IRC[-1]
        
        file_tail=os.path.split(self.file1)[1]        
        filename1=os.path.join(new_dir,file_tail[:-8]+'_React.com')
        filename2=os.path.join(new_dir,file_tail[:-8]+'_Prod.com')
                       
        with open(filename1,'w') as my_file:
                
            my_file.write('%nprocshared=8\n%mem=4GB\n')
            my_file.write('# '+route_section+'\n\n')
            my_file.write('OPT\n\n0 1\n')
            
            for i in range(len(Reactant)):
                my_file.writelines([str(Reactant['Atomic Symbol'][i]),'\t\t',str(Reactant['X'][i]),'\t\t',str(Reactant['Y'][i]),'\t\t',str(Reactant['Z'][i]),'\n'])
                
            my_file.write('\n')  
            
        with open(filename2,'w') as my_file:
                
            my_file.write('%nprocshared=8\n%mem=4GB\n')
            my_file.write('# '+route_section+'\n\n')
            my_file.write('OPT\n\n0 1\n')
            
            for i in range(len(Product)):
                my_file.writelines([str(Product['Atomic Symbol'][i]),'\t\t',str(Product['X'][i]),'\t\t',str(Product['Y'][i]),'\t\t',str(Product['Z'][i]),'\n'])
                
            my_file.write('\n')  

    


class IRC_2F():
    
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
    
    
    def get_Rx_energies_VS_coordinates_dir1(self):
        
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
    
    
    def get_Rx_energies_VS_coordinates_dir2(self):
        
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
    
    
    def get_HF_forces_file1(self):
        
        file=self.file1
        
        nl=self.get_number_of_lines()[0]
        
        ref_word1='Forces (Hartrees/Bohr)'
        
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
                        
                        
                rf_pos=[]
                
                for r in pos2:
                    
                    for s in pos1:
                        
                        if s<r:
                            
                            match=s
                        else:
                            
                            break
                                              
                    rf_pos.append(match)                           
                
                        
                rf_pos.insert(0,pos1[0])
                                
                irc_rf=[]
                
                count_rfm=0
                
                for pos in rf_pos:
                    
                    irc_rf.append([])
                    
                    my_file.seek(pos,0)
                    
                    for i in range(5):
                        
                        my_file.readline()
                        
                    count_atoms=0
                    
                    while True:
                        
                        lc=my_file.readline()
                    
                        if '------' not in lc:
                            
                            irc_rf[count_rfm].append([])
                            
                            irc_rf[count_rfm][count_atoms].append(eval(lc.split()[0]))
                            irc_rf[count_rfm][count_atoms].append(dict_atoms[int(lc.split()[1])])
                            irc_rf[count_rfm][count_atoms].append(eval(lc.split()[2]))
                            irc_rf[count_rfm][count_atoms].append(eval(lc.split()[3]))
                            irc_rf[count_rfm][count_atoms].append(eval(lc.split()[4]))
                            
                            count_atoms+=1
                            
                        else:
                            
                            break
                    
                    
                    count_rfm+=1
            
        else:
            
            print(file,'not found!')
            
        return irc_rf
    
    
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
    
    
    def get_HF_forces_file2(self):
        
        file=self.file2
        
        nl=self.get_number_of_lines()[1]
        
        ref_word1='Forces (Hartrees/Bohr)'
        
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
                        
                        
                rf_pos=[]
                
                for r in pos2:
                    
                    for s in pos1:
                        
                        if s<r:
                            
                            match=s
                        else:
                            
                            break
                                              
                    rf_pos.append(match)                           
                
                        
                rf_pos.insert(0,pos1[0])
                          
                
                irc_rf=[]
                
                count_rfm=0
                
                for pos in rf_pos:
                    
                    irc_rf.append([])
                    
                    my_file.seek(pos,0)
                    
                    for i in range(5):
                        
                        my_file.readline()
                        
                    
                       
                    count_atoms=0
                    
                    while True:
                        
                        lc=my_file.readline()
                    
                        if '------' not in lc:
                            
                            irc_rf[count_rfm].append([])
                            
                            irc_rf[count_rfm][count_atoms].append(eval(lc.split()[0]))
                            irc_rf[count_rfm][count_atoms].append(dict_atoms[int(lc.split()[1])])
                            irc_rf[count_rfm][count_atoms].append(eval(lc.split()[2]))
                            irc_rf[count_rfm][count_atoms].append(eval(lc.split()[3]))
                            irc_rf[count_rfm][count_atoms].append(eval(lc.split()[4]))
                            
                            count_atoms+=1
                            
                        else:
                            
                            break
                    
                    
                    count_rfm+=1
                
            
        else:
            
            print(file,'not found!')
            
            
        
        return irc_rf
        
    
    def construct_geoms_string(self):
        
        geom1=self.get_geometries_file1()
        geom2=self.get_geometries_file2()
        
        dir1=self.get_Rx_energies_VS_coordinates_dir1()
        dir2=self.get_Rx_energies_VS_coordinates_dir2()
        
        gap1=abs(dir1[-1][1]-dir1[0][1])
        gap2=abs(dir2[-1][1]-dir2[0][1])
        
        if gap2>gap1:
            
            geom_dir1=[]
            geom_dir2=[]
            
            if dir1[0][0]>=0 and dir1[1][0]>0:
                
                geom_dir1=self.reverse_list(geom1)
                
            else:
                
                geom_dir1=geom1           
            
            
            if dir2[0][0]<0 and dir2[1][0]<0:
                
                geom_dir2=self.reverse_list(geom2)

            else:
                
                geom_dir2=geom2
                
                
            
            rev=geom_dir1
            forw=geom_dir2
            
            
        else:
            geom_dir1=[]
            geom_dir2=[]
            
            if dir2[0][0]>=0 and dir2[1][0]>0:

                geom_dir2=self.reverse_list(geom2)
                
            else:
                
                geom_dir2=geom2                  
              
            
            if dir1[0][0]<0 and dir1[1][0]<0:
                
                geom_dir1=self.reverse_list(geom1)
                
            else:
                
                geom_dir1=geom1
            
            rev=geom_dir2
            forw=geom_dir1
            
        irc=rev+forw[1:]

        
        return irc

    
    def construct_HF_forces_string(self):
        
        rf1=self.get_HF_forces_file1()
        rf2=self.get_HF_forces_file2()
        
        dir1=self.get_Rx_energies_VS_coordinates_dir1()
        dir2=self.get_Rx_energies_VS_coordinates_dir2()
        
        gap1=abs(dir1[-1][1]-dir1[0][1])
        gap2=abs(dir2[-1][1]-dir2[0][1])
        
        if gap2>gap1:
            
            rf1_dir=[]
            rf2_dir=[]
            
            if dir1[0][0]>=0 and dir1[1][0]>0:
                
                rf1_dir=self.reverse_list(rf1)
                
            else:
                
                rf1_dir=rf1             
            
            
            if dir2[0][0]<0 and dir2[1][0]<0:
                
                rf2_dir=self.reverse_list(rf2)

            else:
                
                rf2_dir=rf2
                
                
            
            rev=rf1_dir
            forw=rf2_dir
            
            
        else:
            rf1_dir=[]
            rf2_dir=[]
            
            if dir2[0][0]>=0 and dir2[1][0]>0:

                rf2_dir=self.reverse_list(rf2)
                
            else:
                
                rf2_dir=rf2                    
              
            
            if dir1[0][0]<0 and dir1[1][0]<0:
                
                rf1_dir=self.reverse_list(rf1)
                
            else:
                
                rf1_dir=rf1
            
            rev=rf2_dir
            forw=rf1_dir
                
        irc=rev+forw[1:]
        
        return irc
                
        
    
    def Df_geoms_string(self):
        
        import pandas as pd
        
        geoms=self.construct_geoms_string()
        
        result=[]
        
        for geom in geoms:
            
            data=pd.DataFrame(geom,columns=['Atomic Symbol','X','Y','Z'])
            
            result.append(data)
    
        return result

    def Df_force_matrix_string(self):
        
        import pandas as pd        
        
        fms=self.construct_HF_forces_string()
        
        result=[]
        
        for fm in fms:
            
            
            data=pd.DataFrame(fm,columns=['Atomic Index','Atomic Symbol','Fx','Fy','Fz'])
            
            result.append(data)
    
        return result           
                      
    
    def reverse_list(self,my_list):
        
        new_list=[]
        
        for i in range(len(my_list)):
            
            new_list.append(my_list[len(my_list)-1-i])
            
        return new_list
        
    def get_Rx_energies_VS_coordinates(self):
        
        dir1=self.get_Rx_energies_VS_coordinates_dir1()
        dir2=self.get_Rx_energies_VS_coordinates_dir2()
        
                
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
                        dir2_new[count].append(dir1[r][1])
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
    
    def direction(self):
        
        en=self.get_Energy()
        
        if abs(en[-1])>abs(en[0]):
            
            result='Normal'
        else:
            result='Reverse'
            
        return result
    
    def get_Ea(self):
        
        En=self.get_Energy()
        
        Ea=En[En.index(max(En))]-En[0]
        
        return Ea

    def get_mid_act_Rx_coord(self):
        
        Energy=self.get_Energy()
        
        Rx=self.get_Rx_coords()
        
        Ea=self.get_Ea()
        
        mid_Rx=Rx[0]
        
        for i in Rx[1:]:
            
            lower_limit=0.50*Ea
            
            diff=Ea=Energy[Rx.index(0)]-Energy[Rx.index(i)]
            
            if diff < Energy[Rx.index(i)]:
                
                mid_Rx=i
                
                break
                
        return mid_Rx, Ea, Rx.index(mid_Rx)
    
    
    def get_mid_act_Rx_coord_2(self):
        
        Energy=self.get_Energy()
        
        Rx=self.get_Rx_coords()
        
        Ea=-(Energy[Rx.index(Rx[-1])]-Energy[Rx.index(0)])
        
        
        mid_Rx=Rx[0]
        
        for i in Rx[Rx.index(0):]:
            
            upper_limit=0.50*Ea
            
            diff=-(Energy[Rx.index(Rx[-1])]-Energy[Rx.index(i)])
            
            
            if diff<upper_limit:
                
                mid_Rx=i
                
                break               
                
                
        return mid_Rx, Ea, Rx.index(mid_Rx)
    
    
                
    def plot_E_VS_Rx_coords(self):
        
        import matplotlib.pyplot as plt
        
        Rc=self.get_Rx_coords()

        Energy=self.get_Energy()
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]

        plt.figure(figsize=(12,8))
        
        plt.ylim(min(Energy)-6,max(Energy)+6)
        
        plt.axvline(alpha_limit,0,color='black')
        
        plt.axvline(gamma_limit,0,color='black')
        
        plt.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plt.axhline(0)
        
        plt.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plt.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plt.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plt.scatter(Rc,Energy,c='blue',linewidths=0.5)
        
        plt.ylabel('Potential Energy, E(ξ)',fontsize=20)
        
        plt.xlabel('Reaction coordinate ξ',fontsize=20)
        
        plt.text((Rc[-1]+gamma_limit+1)/2-2.8,max(Energy)+3,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow(Rc[-1]-gamma_limit,max(Energy)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow(Rc[-1]-gamma_limit,max(Energy)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        plt.text((Rc[0]+alpha_limit)/2-3,max(Energy)+3,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(Energy)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(Energy)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.text((alpha_limit+gamma_limit)/2 - 1.5,max(Energy)+3,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(Energy)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(Energy)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.text(alpha_limit-0.5,min(Energy)-3.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        plt.text(0.4,min(Energy)-3.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plt.text(gamma_limit+0.5,min(Energy)-3.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
        
        
        plt.show()
              

        
    
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
        
        return force,Rc
    
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
                
        for j in range(1,len(force)):
                
            if force[j]>test:
                
                global_max=Rc[j]
                
                test=force[j]
                
        return global_min,global_max
    
    
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
        
    
    
    def plot_Rx_force(self):
        
        import matplotlib.pyplot as plt
        
        Rc=self.get_Rx_coords()
            
        forces=self.get_Rx_force()[0]
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]
        
        plt.figure(figsize=(12,8))
        
        plt.scatter(Rc,forces,c='blue')
        
        plt.ylim(min(forces)-3,max(forces)+3)
        
        plt.axvline(alpha_limit,0,color='black')
        
        plt.axvline(gamma_limit,0,color='black')
        
        #plt.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plt.axhline(0)
        
        plt.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plt.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plt.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plt.ylabel('Reaction force,F(ξ)',fontsize=20)
        
        plt.xlabel('Reaction coordinate (ξ)',fontsize=20)
        
        plt.text((Rc[-1]+gamma_limit+1)/2-2.8,max(forces)+2,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow(Rc[-1]-gamma_limit,max(forces)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow(Rc[-1]-gamma_limit,max(forces)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        plt.text((Rc[0]+alpha_limit)/2-3,max(forces)+2,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(forces)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(forces)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.text((alpha_limit+gamma_limit)/2 - 1.5,max(forces)+2,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(forces)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(forces)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.text(alpha_limit-0.5,min(forces)-2.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        #plt.text(0.4,min(forces)-2.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plt.text(gamma_limit+0.5,min(forces)-2.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
        #plt.savefig('fig.png')
        
        plt.show()
        
    
    
    def plot_Rx_force_constant(self):
        
        import matplotlib.pyplot as plt
        
        Rc=self.get_Rx_coords()
        
        force_constant=self.get_Rx_force_constant()[0]
        
        plt.figure(figsize=(12,8),dpi=600)
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]
        
        plt.ylim(min(force_constant)-3,max(force_constant)+3)
        
        plt.axvline(alpha_limit,0,color='black')
        
        plt.axvline(gamma_limit,0,color='black')
        
        #plt.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plt.axhline(0)
        
        plt.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plt.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plt.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plt.scatter(Rc,force_constant,c='blue')      
        
        plt.ylabel('Reactant force constant,κ(ξ)',fontsize=20)
        
        plt.xlabel('Reaction Coordinates (ξ)',fontsize=20)
        
        plt.text((Rc[-1]+gamma_limit+1)/2-2.8,max(force_constant)+2,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow(Rc[-1]-gamma_limit,max(force_constant)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow(Rc[-1]-gamma_limit,max(force_constant)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        plt.text((Rc[0]+alpha_limit)/2-3,max(force_constant)+2,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(force_constant)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((Rc[0]+alpha_limit)/2,max(force_constant)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.text((alpha_limit+gamma_limit)/2 - 1.5,max(force_constant)+2,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(force_constant)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.arrow((alpha_limit+gamma_limit)/2,max(force_constant)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plt.text(alpha_limit-0.5,min(force_constant)-2.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        #plt.text(0.4,min(force_constant)-2.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plt.text(gamma_limit+0.5,min(force_constant)-2.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
        
        plt.show()
    
    def plot_all(self):
        
        import matplotlib.pyplot as plt
        
        Rc=self.get_Rx_coords()
        
        Energy=self.get_Energy()
        
        forces=self.get_Rx_force()[0]
        
        force_constant=self.get_Rx_force_constant()[0]
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]
        
        plt.figure(figsize=(32,45),dpi=(200))
        
        plot1=plt.subplot2grid((3,2),(0,0),colspan=2)
        
        plot2=plt.subplot2grid((3,2),(1,0),colspan=2)
        
        plot3=plt.subplot2grid((3,2),(2,0),colspan=2)
                
        #Plot 1
        
        plot1.set_ylim(min(Energy)-6,max(Energy)+6)
        
        plot1.axvline(alpha_limit,0,color='black')
        
        plot1.axvline(gamma_limit,0,color='black')
        
        plot1.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plot1.axhline(0)
        
        plot1.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plot1.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plot1.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plot1.scatter(Rc,Energy,c='blue',linewidths=0.5)
        
        plot1.set_ylabel('Potential Energy, E(ξ)',fontsize=20)
        
        plot1.set_xlabel('Reaction coordinate ξ',fontsize=20)
        
        plot1.text(gamma_limit+1,max(Energy)+3,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot1.arrow(Rc[-1]-gamma_limit,max(Energy)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot1.arrow(Rc[-1]-gamma_limit,max(Energy)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        plot1.text(Rc[0]+1,max(Energy)+3,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot1.arrow((Rc[0]+alpha_limit)/2,max(Energy)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot1.arrow((Rc[0]+alpha_limit)/2,max(Energy)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot1.text((alpha_limit+gamma_limit)/2 - 1.5,max(Energy)+3,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot1.arrow((alpha_limit+gamma_limit)/2,max(Energy)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot1.arrow((alpha_limit+gamma_limit)/2,max(Energy)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot1.text(alpha_limit-0.5,min(Energy)-3.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        plot1.text(0.4,min(Energy)-3.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plot1.text(gamma_limit+0.5,min(Energy)-3.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
        
        #Plot2
        
        plot2.set_ylim(min(forces)-3,max(forces)+3)
        
        plot2.axvline(alpha_limit,0,color='black')
        
        plot2.axvline(gamma_limit,0,color='black')
        
        plot2.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plot2.axhline(0)
        
        plot2.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plot2.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plot2.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plot2.scatter(Rc,forces,color='blue')
        
        plot2.set_ylabel('Reaction force,F(ξ)',fontsize=20)
        
        plot2.set_xlabel('Reaction coordinate (ξ)',fontsize=20)
        
        plot2.text(gamma_limit+1,max(forces)+2,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot2.arrow(Rc[-1]-gamma_limit,max(forces)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot2.arrow(Rc[-1]-gamma_limit,max(forces)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        plot2.text(Rc[0]+1,max(forces)+2,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot2.arrow((Rc[0]+alpha_limit)/2,max(forces)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot2.arrow((Rc[0]+alpha_limit)/2,max(forces)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot2.text((alpha_limit+gamma_limit)/2 - 1.5,max(forces)+2,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot2.arrow((alpha_limit+gamma_limit)/2,max(forces)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot2.arrow((alpha_limit+gamma_limit)/2,max(forces)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot2.text(alpha_limit-0.5,min(forces)-2.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        plot2.text(0.4,min(forces)-2.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plot2.text(gamma_limit+0.5,min(forces)-2.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
        
        #Plot 3
        
        plot3.set_ylim(min(force_constant)-3,max(force_constant)+3)
        
        plot3.axvline(alpha_limit,0,color='black')
        
        plot3.axvline(gamma_limit,0,color='black')
        
        plot3.axvline(0,ymax=0.9,color='red',linestyle='--',alpha=0.4)
        
        plot3.axhline(0)
        
        plot3.axvspan(Rc[0],alpha_limit,color='#27c9ae',alpha=0.3)
        
        plot3.axvspan(alpha_limit,gamma_limit,color='#2788c9',alpha=0.3)
        
        plot3.axvspan(gamma_limit,Rc[-1],color='#27c9ae',alpha=0.3)
        
        plot3.scatter(Rc,force_constant,c='blue')      
        
        plot3.set_ylabel('Reactant force constant,κ(ξ)',fontsize=20)
        
        plot3.set_xlabel('Reaction Coordinates (ξ)',fontsize=20)
        
        plot3.text(gamma_limit+1,max(force_constant)+2,'Relaxation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot3.arrow(Rc[-1]-gamma_limit,max(force_constant)+1.5,(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot3.arrow(Rc[-1]-gamma_limit,max(force_constant)+1.5,-(Rc[-1]-gamma_limit-0.5)/2,0,head_length = 0.15, head_width = 0.15, length_includes_head = True)
        
        plot3.text(Rc[0]+1,max(force_constant)+2,'Preparation Phase',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot3.arrow((Rc[0]+alpha_limit)/2,max(force_constant)+1.5,-(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot3.arrow((Rc[0]+alpha_limit)/2,max(force_constant)+1.5,(alpha_limit-0.5-Rc[0])/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot3.text((alpha_limit+gamma_limit)/2 - 1.5,max(force_constant)+2,'TS region',fontsize=14,color='blue',**{'fontname':'Arial'})
        
        plot3.arrow((alpha_limit+gamma_limit)/2,max(force_constant)+1.5,-(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot3.arrow((alpha_limit+gamma_limit)/2,max(force_constant)+1.5,(gamma_limit-alpha_limit-0.5)/2,0,head_length = 0.15, head_width = 0.1, length_includes_head = True)
        
        plot3.text(alpha_limit-0.5,min(force_constant)-2.5,'α',fontsize=14,color='black',**{'fontname':'Arial'})
        
        plot3.text(0.4,min(force_constant)-2.5,'β',fontsize=14,color='black',**{'fontname':'Arial'})      
        
        plot3.text(gamma_limit+0.5,min(force_constant)-2.5,'γ',fontsize=14,color='black',**{'fontname':'Arial'})
        
                
    def get_TS_region_extend(self):
        
        alpha_limit=self.get_min_max_Rx_force()[0]
        
        gamma_limit=self.get_min_max_Rx_force()[1]
        
        return round(abs(gamma_limit-alpha_limit),2)
        
        
    def generate_irc_geom_spcalc_com_files(self,new_dir,level_of_theory):
        
        irc_path=self.construct_geoms_string()
        
        file_tail=os.path.split(self.file1)[1]
        
        count=1
        
        for geom in irc_path:
            
            filename=os.path.join(new_dir,file_tail[:-8]+'_pt_'+str(count)+'.com')
            
            with open(filename,'w') as my_file:
                
                my_file.write('%nprocshared=8\n%mem=4GB\n')
                my_file.write('# '+level_of_theory+'pop=nbo \n\n')
                my_file.write('SP\n\n0 1\n')
                
                for line in geom:
                    
                    my_file.writelines([str(line[0]),'\t\t',str(line[1]),'\t\t',str(line[2]),'\t\t',str(line[3]),'\n'])

                my_file.write('\n')  

            count+=1   
            

    def decompose_activation_energy(self):

        ALPHA=self.get_min_max_Rx_force()[0]
        BETA=0
        GAMMA=self.get_min_max_Rx_force()[1]
        Rc=self.get_Rx_coords()
        Energy=self.get_Energy()
        Alpha_index=Rc.index(ALPHA)
        Beta_index=Rc.index(BETA)
        E_TS=Energy[Beta_index]
        E_Reactant=Energy[0]
        E_Alpha=Energy[Alpha_index]
        
        return round(E_TS-E_Reactant,2),round(E_Alpha-E_Reactant,2),round(E_TS-E_Alpha,2)          

    
    def get_geoms_Dataframes(self):
        
        import pandas as pd
        
        irc_path=self.construct_geoms_string()
        IRC=[]
        for geom in irc_path:
            
            data=pd.DataFrame(geom,columns=['Atomic Symbol','X','Y','Z'])
            data['Atomic Label']=[x for x in range(1,len(data)+1)]
            IRC.append(data)
            
        return IRC 
    
    def get_reduced_geoms_Dataframes(self):  # These are H depleted geometries
    
        IRC=self.get_geoms_Dataframes()
        new_IRC=[]
        
        for geom in IRC:
            
            H_depleted_geom=geom[geom['Atomic Symbol']!='H']
            
            new_IRC.append(H_depleted_geom)
            
        return new_IRC
                
    
    def compute_distance(self,atom1,atom2):
        
        import math
        import numpy as np
        
        coords_zipper=zip(atom1,atom2)
        vect=[(r-s) for r,s in coords_zipper]        
        dist=np.linalg.norm(vect)
        
        return dist
    
    
    
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
                
                
    def get_unique_input(self,output_dir,job_details):
        
        file=self.filename
            
        TS=self.Df_ts_geometry()        
            
        file_tail=os.path.split(file)[1]
            
        new_filename=os.path.join(output_dir,'IRC'+file_tail[2:-3]+'gjf')
            
        with open(new_filename,'w') as myfile:
                
            myfile.write('%nprocshared=8\n%mem=4GB\n')
            
            myfile.write(job_details)
                
            myfile.write('\n\nIRC\n\n0 1\n')
                
            if self.check_convergence()=='Successful':
                    
                for index in TS.index:
                        
                    myfile.writelines([str(dict_atoms[TS['Atomic Number'][index]]),'\t\t',str(self.to_decimal(TS['X'][index])),'\t\t',str(self.to_decimal(TS['Y'][index])),'\t\t',str(self.to_decimal(TS['Z'][index])),'\n'])
                        
                myfile.write('\n')
                    
            else:
                    
                for index in TS.index:
                        
                    myfile.writelines([str([TS['Atomic Symbol'][index]]),'\t\t',str(self.to_decimal(TS['X'][index])),'t\t',str(self.to_decimal(TS['Y'][index])),'\)t\t',str(self.to_decimal(TS['Z'][index])),'\n'])
                        
                myfile.write('\n')
    
    
    
    
    
    
    