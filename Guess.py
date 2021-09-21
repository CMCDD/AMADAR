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
import pickle



dict_atoms={
    
    1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C',
    
    7: 'N', 8:'O', 9:'F', 10:'Ne', 11:'Na', 12: 'Mg',
    
    13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl', 18: 'Ar',
    
    19: 'K', 20: 'Ca', 21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',
    
    27:'Co',28:'Ni',29:'Cu',30:'Zn',31:'Ga',32:'Ge',33:'As',34:'Se',
    
    35:'Br',36:'Kr',37:'Rb',38:'Sr',39:'Y',40:'Zr',41:'Nb',42:'Mo',43:'Tc',
    
    44:'Ru',45:'Rh',46:'Pd',47:'Ag',48:'Cd',49:'In',50:'Sn',51:'Sb',52:'Te',
    
    53:'I',54:'Xe'}


class RCO_guess_ts():   #RCO_pguess_ts()
    
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
            
    def degree_of_freedom(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Deg. of freedom'
        
        if os.path.isfile(file):
            
             with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                        
                        degr=int(lc.split()[3])          
            
            
        else:
        
            print('File Not Found')
            
        return degr
    
    
    def number_of_distinct_conf(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Optimization completed'
        
        count=0
        
        if os.path.isfile(file):
            
             with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                        
                        count+=1
                        
        return count
    
    
    def locate_confs(self):
        
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
        
        list_scf_done_lines_pos=[]
        
        list_stand_orient_lines_pos=[]
        
        list_opt_param_lines_pos=[]
        
        for pos3 in pos3_opt_compl:
            
            scf_macth=0
            
            stand_or_match=0
            
            for pos1 in pos1_scf:
                
                if pos1<pos3:
                    
                    scf_match=pos1
                    
            list_scf_done_lines_pos.append(scf_match-2*len('SCF Done'))
            
            for pos2 in pos2_stand_or:
                
                if pos2<pos3:
                    
                    stand_or_match=pos2
                    
                
            list_stand_orient_lines_pos.append(stand_or_match-len('Standard orientation'))
            
            for pos4 in pos4_opt_param:
                
                if pos4>pos3 and (pos4-pos3)<250:
                    
                    opt_param_match=pos4
                
            list_opt_param_lines_pos.append(opt_param_match-len('!   Optimized Parameters   !'))
            
                    
            
        return list_scf_done_lines_pos,list_stand_orient_lines_pos,list_opt_param_lines_pos
    
    
    def get_opt_steps_energies(self):
        
        #This function locate different stationary point of the PES, the energy of the last conf as well
        
        file=self.filename
        
        scf_pos=self.locate_confs()[0]
        
        list_energies=[]
        
        with open(file,'r') as my_file:
            
            for pos in scf_pos:
            
                my_file.seek(pos)
                
                my_file.readline()
                
                lc=my_file.readline()
                
                
                list_energies.append(float(lc.split()[4]))
                
        
        energy_last_conf=max(list_energies)            
                
                
        return list_energies,energy_last_conf
    
    def plot_opt_steps(self,center_atom_interest_dn1,center_atom_interest_dp1,center_atom_interest_dn2,center_atom_interest_dp2):
        
        import matplotlib.pyplot as plt
        
        import numpy as np
        
        energies=self.get_opt_steps_energies()[0]
        
        y=[(x-min(energies)) for x in energies]
        
        y1=[self.compute_dist_atoms(x,[center_atom_interest_dn1,center_atom_interest_dp1]) for x in range(1,len(energies)+1)]
        
        y2=[self.compute_dist_atoms(x,[center_atom_interest_dn2,center_atom_interest_dp2]) for x in range(1,len(energies)+1)]
        
        scaling_factor=int(max(y1)/max(y))
        
        y1=[x/scaling_factor for x in y1]
        
        y2=[x/scaling_factor for x in y2]
        
        last_conf=max(y)
        
        x=np.arange(1,len(energies)+1)
        
        plt.plot(x,y,c='red',label='relative energy (A.U.)')
        
        plt.plot(x,y1,label='dist1*'+str(scaling_factor)+'(Å)',c='green')
        
        plt.plot(x,y2,label='dist2*'+str(scaling_factor)+'(Å)',c='yellow')
        
        plt.xlabel('Optimization steps')
        
        plt.ylabel('Relative Electronic energy (A.U.)')
        
        plt.legend()
        
        plt.show()
        
    
    def get_confs_geometries(self):
        
        #This function locate geometries of different stationary points 
        
        file=self.filename
        
        stand_or_pos=self.locate_confs()[1]
        
        list_conf_geometries=[]
        
        with open(file,'r') as my_file:
            
            count_pos=0
            
            for pos in stand_or_pos:
                
                list_conf_geometries.append([])
            
                my_file.seek(pos)
                
                my_file.readline()
                
                my_file.readline()
                
                my_file.readline()
                
                my_file.readline()
                
                my_file.readline()
                
                my_file.readline()
                
                lc=my_file.readline()
                
                line_count=0
                
                while '----------' not in lc:
                    
                    list_conf_geometries[count_pos].append([])
                    
                    for i in range(len(lc.split())):
                    
                        list_conf_geometries[count_pos][line_count].append(eval(lc.split()[i]))
                        
                    line_count+=1
                    
                    lc=my_file.readline()
                
                count_pos+=1
                
                
        return list_conf_geometries
    
    
    def Df_conf_geometries(self,index=16):
        
        import pandas as pd
        
        confs=self.get_confs_geometries()
        
        list_data=[]
        
        for conf in confs:
            
            data=pd.DataFrame(conf,columns=['Center Number','Atomic Number','Atomic Type','X','Y','Z'])
            
            list_data.append(data)
            
        
        return list_data[index-1]
    
    
    def conf_opt_bond_distances(self,index=13):
        
        file=self.filename
        
        pos4_opt_param=self.locate_confs()[2]
        
        
        list_conf_opt_param=[]
        
        with open(file,'r') as my_file:
            
            count_pos=0
            
            for pos in pos4_opt_param:
                
                list_conf_opt_param.append([])
                
                my_file.seek(pos,0)
                
                my_file.readline()
                
                my_file.readline()
                 
                my_file.readline()
                
                my_file.readline()
                
                my_file.readline()
                
                my_file.readline()
                
                count_bond=0
                
                while True:
                    
                
                    lc=my_file.readline()
                
                    if 'R' in lc.split()[1]:
                        
                        list_conf_opt_param[count_pos].append([])
                        
                        bond_atoms=(lc.split()[2])[2:-1]
                    
                        begin_atom,end_atom=bond_atoms.split(',')[0],bond_atoms.split(',')[1]
                        
                        bond_length=lc.split()[3]
                        
                        list_conf_opt_param[count_pos][count_bond].append(eval(begin_atom))
                        
                        list_conf_opt_param[count_pos][count_bond].append(eval(end_atom))
                        
                        list_conf_opt_param[count_pos][count_bond].append(eval(bond_length))
                        
                        count_bond+=1
                        
                        
                        
                    else:
                        
                        break
                    
                
                count_pos+=1
                
        return list_conf_opt_param[index-1]
    
    
    def get_atomic_number(self,conf_index,atomic_center):
        
        conf=self.Df_conf_geometries(conf_index)
        
        atomic_number=conf['Atomic Number'][atomic_center-1]
        
        
        return atomic_number
    
    
               
   
    def compute_dist_atoms(self,conf_index,list_indices):
        
        from math import sqrt
        
        import numpy as np
        
        if len(list_indices)==2:
            
        
            conf=self.Df_conf_geometries(conf_index)
        
            atom1_coords=[conf['X'][(list_indices[0])-1],conf['Y'][(list_indices[0])-1],conf['Z'][(list_indices[0])-1]]
        
            atom2_coords=[conf['X'][(list_indices[1])-1],conf['Y'][(list_indices[1])-1],conf['Z'][(list_indices[1])-1]]
        
            dist=np.sqrt(np.sum(np.square([(x1-x2) for x1,x2 in zip(atom1_coords,atom2_coords)])))
            
        else:
            
            print('Distance can be computed only between two points. \nMake sure you provided a list of 2 items as second argument to \'compute_dist_atoms()')
        
        
        return dist
    
    def compute_angle(self,conf_index,list_indices):
        
        import numpy as np
        
        import math
        
        
        index_atom1,index_atom2,index_atom3=list_indices[0],list_indices[1],list_indices[2]
        
        if len(list_indices)==3:
            
        
            conf=self.Df_conf_geometries(conf_index)
        
            atom1_coords=[conf['X'][(list_indices[0])-1],conf['Y'][(list_indices[0])-1],conf['Z'][(list_indices[0])-1]]
        
            atom2_coords=[conf['X'][(list_indices[1])-1],conf['Y'][(list_indices[1])-1],conf['Z'][(list_indices[1])-1]]
            
            atom3_coords=[conf['X'][(list_indices[2])-1],conf['Y'][(list_indices[2])-1],conf['Z'][(list_indices[2])-1]]
            
            vect12=[(r2-r1) for r2,r1 in zip(atom2_coords,atom1_coords)]
            
            vect23=[(r3-r2) for r3,r2 in zip(atom3_coords,atom2_coords)]
           
            vect13=[(r3-r1) for r3,r1 in zip(atom3_coords,atom1_coords)]
            
            norm_vect12=np.linalg.norm(vect12)
            
            norm_vect23=np.linalg.norm(vect23)
            
            norm_vect13=np.linalg.norm(vect13)
            
            scalar_prod_vect12_vect23=np.sum(np.array(vect12)*np.array(vect23))
            
            cos_theta=(scalar_prod_vect12_vect23)/(norm_vect12*norm_vect23)
            
            raw_angle=math.degrees(np.arccos(cos_theta))
            
            #Making sure we can avoid angles having the same cos values. 
            #When these angles are one < to 90 deg and the other < 90 deg.
            
            if cos_theta>0 and norm_vect13>math.sqrt(np.sum([x*x for x in [norm_vect12,norm_vect23]])):
                
                theta=raw_angle+2*(90-raw_angle)
            
            elif cos_theta>0 and norm_vect13<math.sqrt(np.sum([x*x for x in [norm_vect12,norm_vect23]])):
                
                theta=raw_angle 
                
            elif cos_theta<0 and norm_vect13>math.sqrt(np.sum([x*x for x in [norm_vect12,norm_vect23]])):
                
                theta=-(raw_angle+2*(90-raw_angle))
            
            else:
                
                theta=-raw_angle
            
            

            
        else:
            
            print('Angle can be computed only between three points. \nMake sure you provided a list of 3 items as second argument to \'compute_angle()')
        
        
        
        return theta
    
    def to_decimal(self,x,precision:int=10):
        
        return format(x,f".{precision}f").lstrip().rstrip('0')
        
    def extract_last_conf(self,output_dir,conf_index,job_details='# opt=(calcall,ts,noeigentest,maxstep=5) pm6 scf=xqc \n\n'):
            
            #extract_last_conf_com_file ()
        
        conf=self.Df_conf_geometries(conf_index)
        
        tail=os.path.split(self.filename)[1]
        new_filename=os.path.join(output_dir,tail[1:-3]+'com')
               
        
        with open(new_filename,'w') as my_new_file:
            
            my_new_file.write('%nprocshared=4\n')
            
            my_new_file.write(job_details)
            
            my_new_file.write('\n\nOptimization\n\n')
            
            my_new_file.write('0 1\n')
                   
            for index in conf.index:
                
                my_new_file.writelines([dict_atoms[conf['Atomic Number'][index]],'\t\t',str(self.to_decimal(conf['X'][index],8)).rjust(15),'\t\t',str(self.to_decimal(conf['Y'][index],8)).rjust(15),'\t\t',str(self.to_decimal(conf['Z'][index],8)).rjust(15),'\n'])
            
            my_new_file.write('\n')


class Opt_GS_error_conv():
    
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
    
    def locate_last_geom(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Standard orientation'
        
        with open(file,'r') as myfile:
            
            for i in range(nl):
                
                lc=myfile.readline()
                
                if ref_word in lc:
                    
                    stand_or_pos=myfile.tell()-len(lc)
                    
                    
        return stand_or_pos
    
    def locate_input_geom(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word=' Charge =  0 Multiplicity = 1'
        
        with open(file,'r') as myfile:
            
            for i in range(nl):
                
                lc=myfile.readline()
                
                if ref_word in lc:
                    
                    input_or_pos=myfile.tell()-len(lc)
                    
                    break
                    
                    
        return input_or_pos
    
    def get_last_geometry(self):
        
                
        file=self.filename
        
        stand_or_pos=self.locate_last_geom()
        
        last_geom=[]
        
        
        with open(file,'r') as my_file:
                
            my_file.seek(stand_or_pos)
            
            for i in range(6):
                
                my_file.readline()
                
                            
            lc=my_file.readline()
                
            line_count=0
                
            while '----------' not in lc:
                
                last_geom.append([])
                    
                    
                for i in range(len(lc.split())):
                    
                    last_geom[line_count].append(eval(lc.split()[i]))
                        
                line_count+=1
                    
                lc=my_file.readline()
                
                
                
                
        return last_geom
    
    
    def get_input_geometry(self):
        
        
        file=self.filename
        
        pos=self.locate_input_geom()
        
        input_geom=[]
        
        
        with open(file,'r') as my_file:
                
            my_file.seek(pos)
            
            my_file.readline()
                
                            
            lc=my_file.readline()
                
            line_count=0
                
            while len(lc.split())==4:
                
                input_geom.append([])
                    
                    
                for i in range(len(lc.split())):
                    
                    if i==0:
                        
                        input_geom[line_count].append(lc.split()[i])
                        
                    else:
                    
                        input_geom[line_count].append(eval(lc.split()[i]))
                        
                line_count+=1
                    
                lc=my_file.readline()            
                
                
        return input_geom
    
    
    def to_decimal(self,x,precision:int=10):
        
        return format(x,f".{precision}f").lstrip().rstrip('0')
    
    def Df_last_geometry(self):
        
        import pandas as pd
        
        minimum=self.get_last_geometry()
        
        data=pd.DataFrame(minimum,columns=['Center Number','Atomic Number','Atomic Type','X','Y','Z'])

        
        return data
    
    def Df_input_geometry(self):
        
        import pandas as pd
        
        minimum=self.get_input_geometry()
        
        data=pd.DataFrame(minimum,columns=['Atomic Symbol','X','Y','Z'])

        
        return data
    
            
    def extract_gjf_file_last_geom(self,output_dir,job_details):
        
        last_geom=self.Df_last_geometry()
        file_tail=os.path.split(self.filename)[1]
        filename=os.path.join(output_dir,file_tail[1:-3]+'gjf')
        
        with open(filename,'w') as my_new_file:
            
            my_new_file.write('%nprocshared=4\n')
            
            my_new_file.write(job_details)
            
            my_new_file.write('\n\nOpt Freq\n\n0 1\n')
            
                   
            for index in last_geom.index:
                
                
                my_new_file.writelines([dict_atoms[last_geom['Atomic Number'][index]],'\t\t',self.to_decimal(last_geom['X'][index],8).rjust(15),'\t\t',self.to_decimal(last_geom['Y'][index],8).rjust(15),'\t\t',self.to_decimal(last_geom['Z'][index],8).rjust(15),'\n'])
    
            my_new_file.write('\n')
            
    def extract_gjf_file_input_geom(self,output_dir,job_details):
        
        input_geom=self.Df_input_geometry()
        file_tail=os.path.split(self.filename)[1]
        filename=os.path.join(output_dir,file_tail[:-3]+'gjf')
        
        with open(filename,'w') as my_new_file:
            
            my_new_file.write('%nprocshared=4\n')
            
            my_new_file.write(job_details)
            
            my_new_file.write('\n\nOpt Freq\n\n0 1\n')
            
                   
            for index in input_geom.index:
                
                
                my_new_file.writelines([input_geom['Atomic Symbol'][index],'\t\t',self.to_decimal(input_geom['X'][index],8).rjust(15),'\t\t',self.to_decimal(input_geom['Y'][index],8).rjust(15),'\t\t',self.to_decimal(input_geom['Z'][index],8).rjust(15),'\n'])
    
            my_new_file.write('\n')
        
        
    

            