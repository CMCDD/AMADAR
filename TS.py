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


dict_atoms={
    
    1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C',
    
    7: 'N', 8:'O', 9:'F', 10:'Ne', 11:'Na', 12: 'Mg',
    
    13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl', 18: 'Ar',
    
    19: 'K', 20: 'Ca', 21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',
    
    27:'Co',28:'Ni',29:'Cu',30:'Zn',31:'Ga',32:'Ge',33:'As',34:'Se',
    
    35:'Br',36:'Kr',37:'Rb',38:'Sr',39:'Y',40:'Zr',41:'Nb',42:'Mo',43:'Tc',
    
    44:'Ru',45:'Rh',46:'Pd',47:'Ag',48:'Cd',49:'In',50:'Sn',51:'Sb',52:'Te',
    
    53:'I',54:'Xe'}



class Opt_TS():
    
    
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
        
        result=None
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
            
                count=len(my_file.readlines())
                
            result=count  
        else:
            
            print('File Not Found ')
            
        return result
    
    
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
    
    
    def extract_ts_gjf_file(self,output_dir,job_details='# opt=(calcall,ts,noeigentest) freq b3lyp/6-31g(d) scf=xqc'):
        
        conf=self.Df_ts_geometry()
        
        tail=os.path.split(self.filename)[1]
        new_filename=os.path.join(output_dir,tail[1:-3]+'com')        
        
        with open(new_filename,'w') as my_new_file:
            
            my_new_file.write('%nprocshared=4\n')
            my_new_file.write(job_details)
            
            my_new_file.write('\n\nTS\n\n')
            
            my_new_file.write('0 1\n')
            
            if self.check_convergence()=='Successful' or self.locate_last_reached_geom()[0]!=0:
                   
                for index in conf.index:
                    
                    my_new_file.writelines([dict_atoms[conf['Atomic Number'][index]],'\t\t',str(round(conf['X'][index],9)).rjust(15),'\t\t',str(round(conf['Y'][index],9)).rjust(15),'\t\t',str(round(conf['Z'][index],9)).rjust(15),'\n'])
        
                my_new_file.write('\n')
                
            else:
                
                for index in conf.index:
                    
                    my_new_file.writelines([conf['Atomic Symbol'][index],'\t\t',str(conf['X'][index]).rjust(15),'\t\t',str(conf['Y'][index]).rjust(15),'\t\t',str(conf['Z'][index]).rjust(15),'\n'])
        
                my_new_file.write('\n')
                
                
            
    
    def get_total_electronic_energy(self,unit='A.U.'):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='SCF Done:'
        
        last_step=''
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                        
                        last_step=lc
                
            raw_energy=float(last_step.split()[4])
                
            if unit.upper()=='EV':
                    
                energy=raw_energy*27.2114
                    
            elif unit.upper()=='KCAL/MOL':
                    
                energy=raw_energy*627.503
                    
            else:
                    
                energy=raw_energy
                    
        else:
            
            print('File Not Found')
            
        return energy
    
    
    def get_ZPE(self,unit='A.U.'):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Zero-point correction'
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                        
                        raw_ZPE=float(lc.split()[2])
                    
                        
            if unit.upper()=='EV':
                    
                ZPE=raw_ZPE*27.2114
                    
            elif unit.upper()=='KCAL/MOL':
                    
                ZPE=raw_ZPE*627.503
                    
            else:
                    
                ZPE=raw_ZPE           
            
        else:
            
            print('File Not Found')
            
        
        return ZPE
    
    
    def get_GFE(self,unit='A.U.'):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Sum of electronic and thermal Free Energies'
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                        
                        raw_GFE=float(lc.split()[7])
                        
                    
                        
            if unit.upper()=='EV':
                    
                GFE=raw_GFE*27.2114
                    
            elif unit.upper()=='KCAL/MOL':
                    
                GFE=raw_GFE*627.503
                    
            else:
                    
                GFE=raw_GFE           
            
        else:
            
            print('File Not Found')
            
        
        return GFE
    
    def get_Enthalpy(self,unit='A.U.'):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Sum of electronic and thermal Enthalpies'
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                        
                        raw_Enthalpy=float(lc.split()[6])
                    
                        
            if unit.upper()=='EV':
                    
                Enthalpy=raw_Enthalpy*27.2114
                    
            elif unit.upper()=='KCAL/MOL':
                    
                Enthalpy=raw_Enthalpy*627.503
                    
            else:
                    
                Enthalpy=raw_Enthalpy           
            
        else:
            
            print('File Not Found')
            
        
        return Enthalpy
    
    
    def get_Sum_electron_ZPE(self,unit='A.U.'):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Sum of electronic and zero-point Energies'
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                        
                        raw_SEZPE=float(lc.split()[6])
                        
                                            
            if unit.upper()=='EV':
                    
                SEZPE=raw_SEZPE*27.2114
                    
            elif unit.upper()=='KCAL/MOL':
                    
                SEZPE=raw_SEZPE*627.503
                    
            else:
                    
                SEZPE=raw_SEZPE           
            
        else:
            
            print('File Not Found')
            
        
        return SEZPE
    
    def get_Sum_electron_TE(self,unit='A.U.'):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Sum of electronic and thermal Energies'
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                        
                        raw_SETE=float(lc.split()[6])
                        
                                            
            if unit.upper()=='EV':
                    
                SETE=raw_SETE*27.2114
                    
            elif unit.upper()=='KCAL/MOL':
                    
                SETE=raw_SETE*627.503
                    
            else:
                    
                SETE=raw_SETE           
            
        else:
            
            print('File Not Found')
            
        
        return SETE
    
    def get_vibr_deg_freedom(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Deg. of freedom'
        
        list_freq=[]
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                    
                        deg_freedom=lc.split()[3]
                        
                        break
                        
        return int(deg_freedom)       
        
    
    def get_vib_freq(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Frequencies -- '
        
        vibr_deg_freedom=self.get_vibr_deg_freedom()
        
        
        list_freq=[]
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                        
                        for value in lc.split()[2:]:
                            
                            list_freq.append(float(value))
                            
        
        list_freq=list_freq[:vibr_deg_freedom]
                    
            
        return list_freq
    
    
    def Analyze_vibr_freq(self):
        
        #lET US EXAMINE VIBRATIONAL FREQUENCIES
        
        threshold=-100
            
        imag_freq=[]
            
        pseudo_imag_freq=[]
        
        real_freq=[]
            
        list_freq=self.get_vib_freq()
            
        for vib_freq in list_freq:
                        
            if vib_freq<threshold:
                
                imag_freq.append(vib_freq)
                
            if threshold<=vib_freq<0:
                
                pseudo_imag_freq.append(vib_freq)
                
            else:
                
                real_freq.append(vib_freq)           
                
                        
        return imag_freq,pseudo_imag_freq
    
    
    def compute_dist_2atoms(self,list_indices):
        
        from math import sqrt
        
        import numpy as np
        
        if len(list_indices)==2:
            
        
            ts=self.Df_ts_geometry()
        
            atom1_coords=[ts['X'][(list_indices[0])-1],ts['Y'][(list_indices[0])-1],ts['Z'][(list_indices[0])-1]]
        
            atom2_coords=[ts['X'][(list_indices[1])-1],ts['Y'][(list_indices[1])-1],ts['Z'][(list_indices[1])-1]]
        
            dist=np.sqrt(np.sum(np.square([(x1-x2) for x1,x2 in zip(atom1_coords,atom2_coords)])))
            
        else:
            
            print('Distance can be computed only between two points. \nMake sure you provided a list of 2 items as second argument to \'compute_dist_atoms()')
        
        
        return dist
    
    def compute_angle_3atoms(self,list_indices):
        
        import numpy as np
        
        import math
        
        
        index_atom1,index_atom2,index_atom3=list_indices[0],list_indices[1],list_indices[2]
        
        if len(list_indices)==3:
            
        
            ts=self.Df_ts_geometry()
        
            atom1_coords=[ts['X'][index_atom1-1],ts['Y'][index_atom1-1],ts['Z'][index_atom1-1]]
        
            atom2_coords=[ts['X'][index_atom2-1],ts['Y'][index_atom2-1],ts['Z'][index_atom2-1]]
            
            atom3_coords=[ts['X'][index_atom3-1],ts['Y'][index_atom3-1],ts['Z'][index_atom3-1]]
            
            vect12=[(r2-r1) for r2,r1 in zip(atom2_coords,atom1_coords)]
            
            vect23=[(r3-r2) for r3,r2 in zip(atom3_coords,atom2_coords)]
           
            vect13=[(r3-r1) for r3,r1 in zip(atom3_coords,atom1_coords)]
            
            norm_vect12=np.linalg.norm(vect12)
            
            norm_vect23=np.linalg.norm(vect23)
            
            norm_vect13=np.linalg.norm(vect13)
            
            scalar_prod_vect12_vect23=np.sum(np.array(vect12)*np.array(vect23))
            
            cos_theta=(scalar_prod_vect12_vect23)/(norm_vect12*norm_vect23)
            
            raw_angle=math.degrees(np.arccos(cos_theta))
            
            
            if 0<raw_angle<90 and math.pow(norm_vect13,2)<=(math.pow(norm_vect12,2)+math.pow(norm_vect23,2)):
                
                theta=raw_angle
                
            elif 0<raw_angle<90 and math.pow(norm_vect13,2)>(math.pow(norm_vect12,2)+math.pow(norm_vect23,2)):
                
                theta=raw_angle + 2*(90-raw_angle)
                
            elif 90<raw_angle<180 and math.pow(norm_vect13,2)>(math.pow(norm_vect12,2)+math.pow(norm_vect23,2)):
                
                theta=raw_angle
                
            elif 90<raw_angle<180 and math.pow(norm_vect13,2)<(math.pow(norm_vect12,2)+math.pow(norm_vect23,2)):
                
                theta=raw_angle-2*(raw_angle-90)
                
            elif raw_angle<0 and math.pow(norm_vect13,2)<=(math.pow(norm_vect12,2)+math.pow(norm_vect23,2)):
                
                theta=180+raw_angle
                
            else:
                
                theta=raw_angle
            
        else:
            
            print('Angle can be computed only between three points. \nMake sure you provided a list of 3 items as second argument to \'compute_angle()')
        
        
        
        return theta 
    
    
    def get_geom_param(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='!   Optimized Parameters   !'
        
        bonds={}
        
        angles={}
        
        dihedral={}
        
        opt_param_pos=None
        
        with open(file,'r') as my_file:
            
            for i in range(nl):
                
                lc=my_file.readline()
                
                if ref_word in lc:
                    
                    opt_param_pos=my_file.tell()-len(lc)
                    
                    break
                else:
                    
                    None                   
                    
            if opt_param_pos!=None:
                
                my_file.seek(opt_param_pos,0)
            
                for j in range(5):
                
                    my_file.readline()
                
                while True:
                
                    new_lc=my_file.readline()
                
                    if '--------------' not in new_lc:
                    
                        if 'R' in new_lc.split()[1]:
                        
                            bond=new_lc.split()[2]
                        
                            bound_atoms=bond[2:-1].split(',')
                        
                            bond_length=eval(new_lc.split()[3])
                        
                            atom1_index,atom2_index=int(bound_atoms[0]),int(bound_atoms[1])
                        
                            bonds[(atom1_index,atom2_index)]=bond_length
                       
                        
                        elif 'A' in new_lc.split()[1]:
                        
                            angle=new_lc.split()[2]
                        
                            linked_atoms=angle[2:-1].split(',')
                        
                            value=eval(new_lc.split()[3])
                        
                            atom1_index,atom2_index,atom3_index=int(linked_atoms[0]),int(linked_atoms[1]),int(linked_atoms[2])
                        
                            angles[(atom1_index,atom2_index,atom3_index)]=value
                        
                        
                        
                        elif 'D' in new_lc.split()[1]:
                        
                            dih=new_lc.split()[2]
                        
                            linked_atoms=dih[2:-1].split(',')
                        
                            value=eval(new_lc.split()[3])
                        
                            atom1_index,atom2_index,atom3_index,atom4_index=int(linked_atoms[0]),int(linked_atoms[1]),int(linked_atoms[2]),int(linked_atoms[3])
                            
                            #dih=set([atom1_index,atom2_index,atom3_index,atom4_index])
                            dihedral[(atom1_index,atom2_index,atom3_index,atom4_index)]=value
                            #dihedral[dih]=value
                        
                        
                        else:
                        
                            None
                        
                    else:
                    
                        break
                        
                else:
                    
                    print('Convergence failed')
                
        return bonds,angles,dihedral
    
    
    def get_bonds(self):
        
        import pandas as pd
        
        bonds=self.get_geom_param()[0]
        
        keys=bonds.keys()
        
        values=[bonds[x] for x in keys]
        
        data=pd.DataFrame()
        
        data['Bonds']=keys
        
        data['Length']=values
        
        return data
    
    def get_angles(self):
        
        import pandas as pd
        
        angles=self.get_geom_param()[1]
        
        keys=angles.keys()
        
        values=[angles[x] for x in keys]
        
        data=pd.DataFrame()
        
        data['Angles']=keys
        
        data['Value']=values
        
    def get_dihedral(self):
        
        import pandas as pd
        
        dih=self.get_geom_param()[2]
        
        keys=dih.keys()
        
        values=[dih[x] for x in keys]
        
        data=pd.DataFrame()
        
        data['Dihedral']=keys
        
        data['Values']=values
        
        
    def get_AOI_table(self):
        
        import pandas as pd
        
        file=self.filename
        
        file_tail=os.path.split(file)[1]
        
        if 'Rx' in file_tail:
            
            AOI_table=pd.read_csv('AOI_homo.csv')
            
        elif 'hDA' in file_tail:
            
            AOI_table=pd.read_csv('AOI_hetero.csv')
            
        else:
            
            AOI_table=pd.read_csv('AOI_confl.csv')
            
            
        
        rx_IDs=AOI_table['Reaction_ID']
        
        rx_IDs=rx_IDs.tolist()
        
        return AOI_table,rx_IDs
    
    def get_atoms_of_interest(self):
        
        file=self.filename
        
        file_tail=os.path.split(file)[1]
        
        AOI_table=self.get_AOI_table()[0]
        
        rx_IDs=self.get_AOI_table()[1]
        
        current_Rx_ID=file_tail[3:-4]
        
        current_Rx_index=rx_IDs.index(current_Rx_ID)
        
        atom0,atom1=AOI_table['Atom0'][current_Rx_index],AOI_table['Atom1'][current_Rx_index]
        
        atom2,atom3=AOI_table['Atom2'][current_Rx_index],AOI_table['Atom3'][current_Rx_index]
        
        atom4,atom5=AOI_table['Atom4'][current_Rx_index],AOI_table['Atom5'][current_Rx_index]
                
            
        return atom0,atom1,atom2,atom3,atom4,atom5
    
    
    def identify_diene_dineophile_AOI(self):
        
        file=self.filename
        
        file_tail=os.path.split(file)[1]
        AOI_table=self.get_AOI_table()[0]
        rx_IDs=self.get_AOI_table()[1]
        
        current_Rx_ID=file_tail[3:-4]
        current_Rx_index=rx_IDs.index(current_Rx_ID)
        atom0,atom1=AOI_table['Atom0'][current_Rx_index],AOI_table['Atom1'][current_Rx_index]
        atom2,atom3=AOI_table['Atom2'][current_Rx_index],AOI_table['Atom3'][current_Rx_index]
        atom4,atom5=AOI_table['Atom4'][current_Rx_index],AOI_table['Atom5'][current_Rx_index]
        
        diene={atom0,atom1,atom2,atom5}
        dienophile={atom3,atom4}
        
        return diene,dienophile
    
    def locate_virtual_bonds(self):
        
        import re
        
        pattern=re.compile('(\d+)')
        
        file=self.filename
        nl=self.get_number_of_lines()
        ref_word='Add virtual bond connecting'
        list_virt_bonds=[]
        
        with open(file,'r') as my_file:
            
            for i in range(nl):
                
                lc=my_file.readline()
                
                if ref_word in lc:
                    
                    atom1=lc.split()[5]
                    match1=pattern.search(atom1)
                    atom1_index=match1.group()
                    atom2=lc.split()[7]
                    match2=pattern.search(atom2)
                    atom2_index=match2.group()
                    list_cur_idx=[int(atom1_index),int(atom2_index)]
                    virt_bond=set(sorted(list_cur_idx))
                    list_virt_bonds.append(virt_bond)
                    
        return list_virt_bonds
        
    
    def locate_fragments(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        
        virt_bonds=self.locate_virtual_bonds()
        
        #print(virt_bonds)
        
        ref_word='Initial Parameters'
        
        pos=[]
        
        list_all_bound_atoms=[]
        
        dict_bonding_tree={}
        
        with open(file,'r') as my_file:
            
            for i in range(nl):
                
                lc=my_file.readline()
                
                if ref_word in lc:
                    
                    pos.append(my_file.tell())
                    
            input_geom_pos=pos[0]
            
            my_file.seek(input_geom_pos,0)
            
            for j in range(4):
                
                my_file.readline()
                
            count=0
                
            while True:
                
                new_lc=my_file.readline()
                
                if '-------------------' not in new_lc:
                    
                    geom_par=new_lc.split()[2]
                    
                    if geom_par[0]=='R':
                        
                        list_all_bound_atoms.append([])
                        
                        paired_atoms=geom_par[2:-1].split(',')
                        
                        atom1=paired_atoms[0]
                        
                        atom2=paired_atoms[1]
                        
                                             
                        if set([int(atom1),int(atom2)]) not in virt_bonds:
                            
                            list_all_bound_atoms[count].append(int(atom1))
                            
                            list_all_bound_atoms[count].append(int(atom2))
                            
                        else:
                            
                            None
                        
                        count+=1
                        
                    else:
                        
                        break
                        
                                            
                else:
                    
                    break
                
        new_list_all_bound_atoms=[]
        
        for item in list_all_bound_atoms:
            
            if len(item)==2:
                
                new_list_all_bound_atoms.append(item)
                
            else:
                
                None
        
        #print(new_list_all_bound_atoms)
            
        
        all_atoms=[]
        
        for m in range(len(new_list_all_bound_atoms)):
            
            
            for i in range(2):
                
                all_atoms.append(new_list_all_bound_atoms[m][i])
                
        all_atoms=set(all_atoms)
        
        
        
        for k in range(len(new_list_all_bound_atoms)):
                
            first_atom=new_list_all_bound_atoms[k][0]
                
            second_atom=new_list_all_bound_atoms[k][1]
                
                
            if first_atom in dict_bonding_tree.keys():
                    
                first_atom_old_value=dict_bonding_tree[first_atom]
                    
                first_atom_new_value=first_atom_old_value.union({second_atom})
                    
                dict_bonding_tree[first_atom]=first_atom_new_value
                
            
            else:
                    
                dict_bonding_tree[first_atom]={second_atom}
                dict_bonding_tree[second_atom]={first_atom}
                
                
            if second_atom in dict_bonding_tree.keys():
                    
                second_atom_old_value=dict_bonding_tree[second_atom]
                    
                second_atom_new_value=second_atom_old_value.union({first_atom})
                    
                dict_bonding_tree[second_atom]=second_atom_new_value
                
            else:
                    
                dict_bonding_tree[second_atom]={first_atom}
                    
                
            
                
                    
        new_bonding_tree={}
        
        for key in dict_bonding_tree:
            
            new_value=sorted(list(dict_bonding_tree[key]))
            
            new_bonding_tree[key]=new_value
            
        #print(new_bonding_tree)   
        
        frag1=[1]
        
        frag2=[]
        
        """
        frag1.append(1)
        
        while True:
            
            for r in new_bonding_tree.keys():
                
                #print(r)
                
                if r in frag1:
                    
                    for s in new_bonding_tree[r]:                        
                                           
                        frag1.append(s)
                    
                    #print(frag)
                    
            break
        
        frag1=list(set(sorted(frag1)))
        
        """
        list_keys=new_bonding_tree.keys()
        
        
        count=0
        
        parsed_keys=[]
        
        while count<len(new_bonding_tree.keys()):
            
            keys=[]
            
            for x in frag1:
                
                if x not in parsed_keys:
                    
                    keys.append(x)
                    
            #print(keys)
                    
            if len(keys)>=1:
                
                curr_key=keys[0]
            
                for r in new_bonding_tree[curr_key]:
                
                    frag1.append(r)
                
                    parsed_keys.append(curr_key)
                    
                    #frag1.append(new_bonding_tree[value])
                    
                count+=1
                
            else:
                break
            
        for p in range(1,len(all_atoms)+1):
            
            if p not in frag1:
                
                frag2.append(p)
                
            else:
                
                None
                
        frag1=sorted(list(set(frag1)))
        frag2=sorted(frag2)
               
        return frag1,frag2
    
    def assign_nature_frag(self):
        
        diene_aoi=self.identify_diene_dineophile_AOI()[0]
        dienophile_aoi=self.identify_diene_dineophile_AOI()[1]
        frag1=self.locate_fragments()[0]
        frag2=self.locate_fragments()[1]
        
        if diene_aoi in frag1:
            
            result=['Dn','Dp']
            
        else:
            
            result=['Dp','Dn']
            
        
        return result
            
    
            
    def extract_fragments(self):
        
        TS=self.Df_ts_geometry()
        f1_atoms=self.locate_fragments()[0]
        f2_atoms=self.locate_fragments()[1]
        
        
        if len(f2_atoms)>1:
            
            frag1=[]
            frag2=[]
            
            for i in f1_atoms:
                center_num=TS['Center Number'].tolist()
                index=center_num.index(i)
                atomic_num=TS['Atomic Number'][index]
                x=TS['X'][index]
                y=TS['Y'][index]
                z=TS['Z'][index]
                
                curr_atom=[dict_atoms[atomic_num],x,y,z]
                
                frag1.append(curr_atom)
                
            for j in f2_atoms:
                
                center_num=TS['Center Number'].tolist()
                index=center_num.index(j)
                atomic_num=TS['Atomic Number'][index]
                x=TS['X'][index]
                y=TS['Y'][index]
                z=TS['Z'][index]
                
                curr_atom=[dict_atoms[atomic_num],x,y,z]
                
                frag2.append(curr_atom)
                
            result=(frag1,frag2)
            
        else:
            
            frag1=[]
            
            for i in f1_atoms:
                center_num=TS['Center Number'].tolist()
                index=center_num.index(i)
                atomic_num=TS['Atomic Number'][index]
                x=TS['X'][index]
                y=TS['Y'][index]
                z=TS['Z'][index]
                
                curr_atom=[dict_atoms[atomic_num],x,y,z]
                
                frag1.append(curr_atom)
            
            result=frag1
            
            
        return result
    
    
    def get_TS_NBO_gjf_file(self,output_dir,job_details):
        
        
        TS=self.Df_ts_geometry()
        
        file=os.path.split(self.filename)[1]
        
        filename=os.path.join(output_dir,file[:-3]+'gjf')
        
        with open(filename,'w') as myfile:
            
            myfile.write('%nprocshared=8\n%mem=4GB\n')
            
            myfile.write(job_details)
            
            myfile.write('\n\nSP_NBO\n\n0 1\n')
            
            
            if self.check_convergence()=='Successful' or self.locate_last_reached_geom()[0]!=0:
                   
                for index in TS.index:
                    
                    myfile.writelines([dict_atoms[TS['Atomic Number'][index]],'\t\t',str(round(TS['X'][index],9)).rjust(15),'\t\t',str(round(TS['Y'][index],9)).rjust(15),'\t\t',str(round(TS['Z'][index],9)).rjust(15),'\n'])
        
                myfile.write('\n$nbo bndidx $end\n')
                
            else:
                
                for index in TS.index:
                    
                    myfile.writelines([TS['Atomic Symbol'][index],'\t\t',str(TS['X'][index]).rjust(15),'\t\t',str(TS['Y'][index]).rjust(15),'\t\t',str(TS['Z'][index]).rjust(15),'\n'])
        
                myfile.write('\n$nbo bndidx $end\n')
            
                   
            
    
    def generate_fragments_com_files(self,path_file1,path_file2,output_dir):
        
        file1=path_file1
        file2=path_file2
        dict_atoms1={}
        dict_atoms2={}
        reactant1_mark=set()
        reactant2_mark=set()
        f1=self.extract_fragments()[0]
        f2=self.extract_fragments()[1]
        #print(f1)
        #print(f2)
        dict_frag1={}
        dict_frag2={}
        frag1_mark=set()
        frag2_mark=set()
        
        ref_word='Symbolic Z-matrix:'
        
        with open(file1,'r') as my_file1:
            
            while True:
                
                lc=my_file1.readline()
                
                if ref_word in lc:
                    
                    pos1=my_file1.tell()
                    
                    break
                
            my_file1.seek(pos1,0)
            
            while True:
                
                lc=my_file1.readline()
                
                if len(lc.split())==4:
                    
                    curr_atom=lc.split()
                    
                    if curr_atom[0] in dict_atoms1.keys():
                        
                        dict_atoms1[curr_atom[0]]+=1
                        
                    else:
                        
                        dict_atoms1[curr_atom[0]]=1
                        
                else:
                    break
                
            
            
            for key in dict_atoms1.keys():
                
                new_term=set({key,dict_atoms1[key]})
                
                reactant1_mark=reactant1_mark.union(new_term)
                
                
        with open(file2,'r') as my_file2:
            
            while True:
                
                lc=my_file2.readline()
                
                if ref_word in lc:
                    
                    pos1=my_file2.tell()
                    
                    break
                
            my_file2.seek(pos1,0)
            
            while True:
                
                lc=my_file2.readline()
                
                if len(lc.split())==4:
                    
                    curr_atom=lc.split()
                    
                    if curr_atom[0] in dict_atoms2.keys():
                        
                        dict_atoms2[curr_atom[0]]+=1
                        
                    else:
                        
                        dict_atoms2[curr_atom[0]]=1
                        
                else:
                    break
                
            
            
            for key in dict_atoms2.keys():
                
                new_term=set({key,dict_atoms2[key]})
                
                reactant2_mark=reactant2_mark.union(new_term)
                
        for i in range(len(f1)):
            
            atom=f1[i][0]
            
            if atom in dict_frag1.keys():
                
                dict_frag1[atom]+=1
                
            else:
                
                 dict_frag1[atom]=1
                 
        for key in dict_frag1.keys():
            
            new_term=set({key,dict_frag1[key]})
            
            frag1_mark=frag1_mark.union(new_term)
            
            
        for j in range(len(f2)):
            
            atom=f2[j][0]
            
            if atom in dict_frag2.keys():
                
                dict_frag2[atom]+=1
                
            else:
                
                 dict_frag2[atom]=1
                 
        for key in dict_frag2.keys():
            
            new_term=set({key,dict_frag2[key]})
            
            frag2_mark=frag2_mark.union(new_term)
            
            
        if frag1_mark==reactant1_mark:
            
            nature_frag1=self.assign_nature_frag()[0]
            nature_frag2=self.assign_nature_frag()[1]
            
            reactant1_filename_tail=os.path.split(file1)[1]
            reactant2_filename_tail=os.path.split(file2)[1]
            
            output_filename1=nature_frag1+reactant1_filename_tail[:-3]+'com'
            output_filename2=nature_frag2+reactant2_filename_tail[:-3]+'com'
            
            with open(os.path.join(output_dir,output_filename1),'w') as out_file1, open(os.path.join(output_dir,output_filename2),'w') as out_file2:
                
                out_file1.write('%nprocshared=8\n%mem=4GB\n# b3lyp/6-31G(d) \n\nSingle Point\n\n0 1\n')
                out_file2.write('%nprocshared=8\n%mem=4GB\n# b3lyp/6-31G(d) \n\nSingle Point\n\n0 1\n')
                
                for k in range(len(f1)):
                    
                    cur_atom=f1[k]
                    
                    out_file1.writelines([str(cur_atom[0]),'\t\t',str(cur_atom[1]),'\t\t',str(cur_atom[2]),'\t\t',str(cur_atom[3]),'\n'])
                    
                out_file1.write('\n')
                
                
                for m in range(len(f2)):
                    
                    cur_atom=f2[m]
                    
                    out_file2.writelines([str(cur_atom[0]),'\t\t',str(cur_atom[1]),'\t\t',str(cur_atom[2]),'\t\t',str(cur_atom[3]),'\n'])
                    
                out_file2.write('\n')
                
        else:
            
            nature_frag1=self.assign_nature_frag()[0]
            nature_frag2=self.assign_nature_frag()[1]
            
            reactant1_filename_tail=os.path.split(file1)[1]
            reactant2_filename_tail=os.path.split(file2)[1]
            
            output_filename1=nature_frag1+reactant2_filename_tail[:-3]+'com'
            output_filename2=nature_frag2+reactant1_filename_tail[:-3]+'com'
            
            with open(os.path.join(output_dir,output_filename1),'w') as out_file1, open(os.path.join(output_dir,output_filename2),'w') as out_file2:
                
                out_file1.write('%nprocshared=8\n%mem=4GB\n# b3lyp/6-31G(d) \n\nSingle Point\n\n0 1\n')
                out_file2.write('%nprocshared=8\n%mem=4GB\n# b3lyp/6-31G(d) \n\nSingle Point\n\n0 1\n')
                
                for k in range(len(f1)):
                    
                    cur_atom=f1[k]
                    
                    out_file1.writelines([str(cur_atom[0]),'\t\t',str(cur_atom[1]),'\t\t',str(cur_atom[2]),'\t\t',str(cur_atom[3]),'\n'])
                    
                out_file1.write('\n')
                
                
                for m in range(len(f2)):
                    
                    cur_atom=f2[m]
                    
                    out_file2.writelines([str(cur_atom[0]),'\t\t',str(cur_atom[1]),'\t\t',str(cur_atom[2]),'\t\t',str(cur_atom[3]),'\n'])
                    
                out_file2.write('\n')
                
    def get_number_of_fragments(self):
        
        fragments=self.locate_fragments()
        frag1=fragments[0]
        frag2=fragments[1]
        
        
        if len(frag2)>1:
            result=2
        else:
            result=1
            
        return result
            
    def generate_distorted_reactants_com_files(self,output_dir):
        
        file=self.filename
        file_tail=os.path.split(file)[1]
        
        num_frag=self.get_number_of_fragments()
        
        
        if num_frag==2:
            
            fragments=self.extract_fragments()
            frag1=fragments[0]
            frag2=fragments[1]
            frag_nature=self.assign_nature_frag()
            nat_frag1=frag_nature[0]
            nat_frag2=frag_nature[1]
            output_filename1=file_tail[3:-4]+nat_frag1+'dist.com'
            output_filename2=file_tail[3:-4]+nat_frag2+'dist.com'
            
            with open(os.path.join(output_dir,output_filename1),'w') as out_file1, open(os.path.join(output_dir,output_filename2),'w') as out_file2:
                
                out_file1.write('%nprocshared=8\n%mem=4GB\n# b3lyp/6-31G(d) \n\nSingle Point\n\n0 1\n')
                out_file2.write('%nprocshared=8\n%mem=4GB\n# b3lyp/6-31G(d) \n\nSingle Point\n\n0 1\n')
                
                for k in range(len(frag1)):
                    
                    cur_atom=frag1[k]
                    
                    out_file1.writelines([str(cur_atom[0]),'\t\t',str(cur_atom[1]),'\t\t',str(cur_atom[2]),'\t\t',str(cur_atom[3]),'\n'])
                    
                out_file1.write('\n')
                
                
                for m in range(len(frag2)):
                    
                    cur_atom=frag2[m]
                    
                    out_file2.writelines([str(cur_atom[0]),'\t\t',str(cur_atom[1]),'\t\t',str(cur_atom[2]),'\t\t',str(cur_atom[3]),'\n'])
                    
                out_file2.write('\n')
            
        else:
            
            frag=self.extract_fragments()
            nature='DD'
            output_filename=file_tail[3:-4]+nature+'dist.com'
            
            with open(os.path.join(output_dir,output_filename),'w') as out_file:
                
                out_file.write('%nprocshared=8\n%mem=4GB\n# b3lyp/6-31G(d) \n\nSingle Point\n\n0 1\n')
                
                
                for k in range(len(frag)):
                    
                    cur_atom=frag[k]
                    
                    out_file.writelines([str(cur_atom[0]),'\t\t',str(cur_atom[1]),'\t\t',str(cur_atom[2]),'\t\t',str(cur_atom[3]),'\n'])
                    
                out_file.write('\n')
                
    
    
    def generate_minima_reactants_PM6_com_files(self,output_dir):
        
        file=self.filename
        file_tail=os.path.split(file)[1]
        
        num_frag=self.get_number_of_fragments()
        
        
        if num_frag==2:
            
            fragments=self.extract_fragments()
            frag1=fragments[0]
            frag2=fragments[1]
            frag_nature=self.assign_nature_frag()
            nat_frag1=frag_nature[0]
            nat_frag2=frag_nature[1]
            output_filename1=file_tail[3:-4]+nat_frag1+'m.com'
            output_filename2=file_tail[3:-4]+nat_frag2+'m.com'
            
            with open(os.path.join(output_dir,output_filename1),'w') as out_file1, open(os.path.join(output_dir,output_filename2),'w') as out_file2:
                
                out_file1.write('# opt freq pm6 scf=tight\n\nOpt\n\n0 1\n')
                out_file2.write('# opt freq pm6 scf=tight\n\nOpt\n\n0 1\n')
                
                for k in range(len(frag1)):
                    
                    cur_atom=frag1[k]
                    
                    out_file1.writelines([str(cur_atom[0]),'\t\t',str(cur_atom[1]),'\t\t',str(cur_atom[2]),'\t\t',str(cur_atom[3]),'\n'])
                    
                out_file1.write('\n')
                
                
                for m in range(len(frag2)):
                    
                    cur_atom=frag2[m]
                    
                    out_file2.writelines([str(cur_atom[0]),'\t\t',str(cur_atom[1]),'\t\t',str(cur_atom[2]),'\t\t',str(cur_atom[3]),'\n'])
                    
                out_file2.write('\n')
            
        else:
            
            frag=self.extract_fragments()
            nature='DD'
            output_filename=file_tail[3:-4]+nature+'m.com'
            
            with open(os.path.join(output_dir,output_filename),'w') as out_file:
                
                out_file.write('# opt freq pm6 scf=tight\n\nOpt\n\n0 1\n')
                
                
                for k in range(len(frag)):
                    
                    cur_atom=frag[k]
                    
                    out_file.writelines([str(cur_atom[0]),'\t\t',str(cur_atom[1]),'\t\t',str(cur_atom[2]),'\t\t',str(cur_atom[3]),'\n'])
                    
                out_file.write('\n')
            

            
            
            
            
