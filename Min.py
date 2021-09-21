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


dict_atoms={
    
    1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C',
    
    7: 'N', 8:'O', 9:'F', 10:'Ne', 11:'Na', 12: 'Mg',
    
    13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl', 18: 'Ar',
    
    19: 'K', 20: 'Ca', 21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',
    
    27:'Co',28:'Ni',29:'Cu',30:'Zn',31:'Ga',32:'Ge',33:'As',34:'Se',
    
    35:'Br',36:'Kr',37:'Rb',38:'Sr',39:'Y',40:'Zr',41:'Nb',42:'Mo',43:'Tc',
    
    44:'Ru',45:'Rh',46:'Pd',47:'Ag',48:'Cd',49:'In',50:'Sn',51:'Sb',52:'Te',
    
    53:'I',54:'Xe'}



class Opt_Freq():

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
    
    def get_level_of_theory(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='# opt freq'
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                        
                        slc=lc.split()
                        
                        level_of_theory=slc[3]
                        
                        method=level_of_theory.split('/')[0]
                        
                        basis=level_of_theory.split('/')[1]
                        
                    else:
                        
                        None
                        
        else:
            
            None
            
        
        return [level_of_theory,method,basis]
    
    def check_convergence(self):
        
        file=self.filename
        
        ref_word='Error'
        
        ref_word2='Normal termination'
        
        gs=open(file,'r')
        
        content=gs.readlines()
        
        result=-1
        
        for line in content:
            
                        
            if ref_word in line and ('Error in corrector energy' not in line):
                
                result=0
                
                break
            
            if ref_word2 in line:
                
                result=1
                
                break
            
        return result
    
    def method_of_calculation(self):
        
        level=self.get_level_of_theory()
        
        return level[1]
    
    def basis_set(self):
        
        level=self.get_level_of_theory()
        
        return level[2]
    
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
    
    
    def get_vib_freq(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
        ref_word='Frequencies -- '
        
        list_freq=[]
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word in lc:
                        
                        for value in lc.split()[2:]:
                            
                            list_freq.append(eval(value))
                            
            
            #lET US EXAMINE VIBRATIONAL FREQUENCIES
            
            imag_freq=[]
            
            real_freq=[]
            
            for vib_freq in list_freq:
                
                if vib_freq<0.00:
                    
                    imag_freq.append(vib_freq)
                    
                else:
                    
                    real_freq.append(vib_freq)
                    
            dict_freq_deg={}  #We are creating a dictionary of vibr.frequencies and their degeneracy
            
            for vib_freq in list_freq:
                
                if vib_freq not in dict_freq_deg:
                    
                    dict_freq_deg[vib_freq]=1
                    
                else:
                    
                    dict_freq_deg[vib_freq]+=1
                    
                    
            
        return list_freq,dict_freq_deg,len(imag_freq),len(real_freq)
    
    def vibr_analysis(self):
        
        imag_freq=self.get_vib_freq()[2]
        
        if imag_freq!=0:
           
           result='Not a minimum'
           
        else:
            
            result='Minimum'
            
        return result
    
    def get_lowest_vibr_freq(self):
        
        vibr_freqs=self.get_vib_freq()[0]
        
        return min(vibr_freqs)
            
            


        
    def get_HOMO_LUMO(self):

        file=self.filename

        ref_word1='Alpha  occ. eigenvalues'

        ref_word2='Alpha virt. eigenvalues'

        nl=self.get_number_of_lines()

        OMO=None

        UMO=None


        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word1 in lc:

                        OMO=lc.split()[len(lc.split())-1]
                
                my_file.seek(0,0)
                
                for j in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word2 in lc:
                        

                        UMO=lc.split()[4]

                        break

                        

        else:

            print('File NOt found')

        HOMO=eval(OMO)

        LUMO=eval(UMO)


        return HOMO,LUMO

    
    def get_chemical_potential(self,unit='kcal/mol'):

        HOMO=self.get_HOMO_LUMO()[0]

        LUMO=self.get_HOMO_LUMO()[1]

        result=(HOMO+LUMO)/2

        if unit.upper()=='KCAL/MOL':

            chem_pot=627.509*result

        elif unit.upper()=='EV':

            chem_pot=27.211*result

        elif unit.upper()=='CM-1':

            chem_pot=219474.6*result

        elif unit.upper()=='KJ/MOL':

            chem_pot=2625.5*result

        else:

            None


        return round(chem_pot,4)


    def get_hardness(self,unit='kcal/mol'):

        HOMO=self.get_HOMO_LUMO()[0]

        LUMO=self.get_HOMO_LUMO()[1]

        result=(LUMO-HOMO)

        if unit.upper()=='KCAL/MOL':

            hardness=627.509*result

        elif unit.upper()=='EV':

            hardness=27.211*result

        elif unit.upper()=='CM-1':

            hardness=219474.6*result

        elif unit.upper()=='KJ/MOL':

            hardness=2625.5*result

        else:

            None



        return round(hardness,4)


    def get_electrophilicity(self,unit='kcal/mol'):

        chem_pot=self.get_chemical_potential(unit)

        hardness=self.get_hardness(unit)


        result= (chem_pot*chem_pot)/(2*hardness)


        return result
    
    
    def locate_minima(self):
        
        file=self.filename
        
        nl=self.get_number_of_lines()
        
               
        ref_word1='Standard orientation'
        
        ref_word2='Optimization completed'
        
        pos1_stand_or=[]
        
        pos2_opt_compl=[]
        
       
        
        if os.path.isfile(file):
            
            with open(file,'r') as my_file:
                
                for i in range(nl):
                    
                    lc=my_file.readline()
                    
                    if ref_word1 in lc:
                        
                        end_line_pos=my_file.tell()
                        
                        pos1_stand_or.append(end_line_pos-len(lc))
                        
                        
                    elif ref_word2 in lc:
                        
                        end_line_pos=my_file.tell()
                        
                        pos2_opt_compl.append(end_line_pos-len(lc))
                        
                    
                    else:
                        
                        None
            
            
        else:
            
            print('File Not Found')
            
        #Let us locate different configurations
        
        stand_orient_line_pos=[]
        
        opt_param_line_pos=[]
        
        for pos2 in pos2_opt_compl:
            
                        
            for pos1 in pos1_stand_or:

                
                if pos1<pos2:
                    
                    stand_or_match=pos1
                    
                
            stand_orient_line_pos=stand_or_match-len('Standard orientation')
            
            
            
        return stand_orient_line_pos
    
    
    def get_minimum_geometry(self):
        
        #This function locate geometries of different stationary points 
        
        file=self.filename
        
        stand_or_pos=self.locate_minima()
        
        
        
        minimum_geom=[]
        
        
        with open(file,'r') as my_file:
                
            my_file.seek(stand_or_pos)
            
            for i in range(6):
                
                my_file.readline()
                
                            
            lc=my_file.readline()
                
            line_count=0
                
            while '----------' not in lc:
                
                minimum_geom.append([])
                    
                    
                for i in range(len(lc.split())):
                    
                    minimum_geom[line_count].append(eval(lc.split()[i]))
                        
                line_count+=1
                    
                lc=my_file.readline()
                
                
                
                
        return minimum_geom
    
    def to_decimal(self,x,precision:int=10):
        
        return format(x,f".{precision}f").lstrip().rstrip('0')
    
    def Df_minimum_geometry(self):
        
        import pandas as pd
        
        minimum=self.get_minimum_geometry()
        
        data=pd.DataFrame(minimum,columns=['Center Number','Atomic Number','Atomic Type','X','Y','Z'])

            
        
        return data
    
    
    def extract_xyz_file_reactants(self,output_dir):
        
        minimum=self.Df_minimum_geometry()
        size=len(minimum)
        file_tail=os.path.split(self.filename)[1]
        filename=os.path.join(output_dir,file_tail[:-3]+'xyz')
        
        with open(filename,'w') as my_new_file:
            
            my_new_file.write(str(size))
            
            my_new_file.write('XYZ')
            
                   
            for index in minimum.index:
                
                my_new_file.writelines([dict_atoms[minimum['Atomic Number'][index]],'\t\t',str(round(minimum['X'][index],9)).rjust(15),'\t\t',str(round(minimum['Y'][index],9)).rjust(15),'\t\t',str(round(minimum['Z'][index],9)).rjust(15),'\n'])
    
            #my_new_file.write('\n')
            
    def extract_gjf_file(self,output_dir,job_details):
        
        minimum=self.Df_minimum_geometry()
        size=len(minimum)
        file_tail=os.path.split(self.filename)[1]
        filename=os.path.join(output_dir,file_tail[:-3]+'com')
        
        with open(filename,'w') as my_new_file:
            
            my_new_file.write('%nprocshared=4\n')
            
            my_new_file.write(job_details)
            
            my_new_file.write('\n\nOpt Freq\n\n0 1\n')
            
                   
            for index in minimum.index:
                
                
                my_new_file.writelines([dict_atoms[minimum['Atomic Number'][index]],'\t\t',self.to_decimal(minimum['X'][index],8).rjust(15),'\t\t',self.to_decimal(minimum['Y'][index],8).rjust(15),'\t\t',self.to_decimal(minimum['Z'][index],8).rjust(15),'\n'])
    
            my_new_file.write('\n')