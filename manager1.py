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
from multiprocessing import Pool
from Gaussian import RunGaussian
from Min import Opt_Freq
import shutil
import time





def Optimize_R(my_path_1,levels,cwd): 
    
    try:
        path_level2=os.path.join(my_path_1,'REFINED')
        errors_path=os.path.join(my_path_1,'ERROR_FILES')
        make_path(path_level2)
        make_path(errors_path)
        errors_path2=os.path.join(path_level2,'ERROR_FILES')
        make_path(errors_path2)
        list_input_files_level1=[]
        list_input_files_level2=[]
        list_output_files_level1=[]
        list_output_files_level2=[]
        
        if os.path.isdir(my_path_1):
            for file in os.listdir(my_path_1):
                if file[-3:]=='com':
                    filepath=os.path.join(my_path_1,file)
                    list_input_files_level1.append(filepath)
        
        print('\nGaussian called for calculation\nSystems: Reactant(s)\nJob Type: OPT FREQ\nLEVEL: ',str(levels[0]).upper(),'\n')
        pool1=Pool(processes=4)
        result=pool1.map(RunGaussian,list_input_files_level1)
        pool1.close()
        pool1.join()   # The join() method block the execution until all the jobs fed into pool1 are terminated!!!      
        print('>>>> Gaussian Execution Terminated')
        
        #Identifying output files that converged and creating new inputs
        normal_conv=[]
        error_conv=[]
        for file in os.listdir(my_path_1):
            if file[-3:]=='log':
                list_output_files_level1.append(file)
                filename=os.path.join(my_path_1,file)
                INSTANCE=Opt_Freq(filename)                
                if INSTANCE.check_convergence()==1:
                    normal_conv.append(filename)
                    job_details='# opt freq '+str(levels[1])+' scf=(xqc, maxcycle=1024)'
                    INSTANCE.extract_gjf_file(path_level2,job_details)
                    shutil.copy(filename,path_level2)
                else:
                    error_conv.append(filename)
                    shutil.copy(filename,errors_path)
        
        print('>>>> Reactants optimized :',str(len(list_output_files_level1)))       
        
        if len(normal_conv)>=1:
            print('>>>> Normally converged :',len(normal_conv))   
            for filename in normal_conv:
                print('    >>>>',filename, '\tSize(kB) :',int((os.stat(filename).st_size)/1024))
                time.sleep(1)
            print('\n    >>>> Structure(s) to be refined. Last geometries saved in :'+str(path_level2)+'>>>>')
        
        if len(error_conv)>=1:
            print('>>>> Error termination or job killed :',len(error_conv))
            for filename in error_conv:
                error_type=get_error_type(filename)
                print('    >>>>',filename, '\tSize(kB) :',int((os.stat(filename).st_size)/1024),' Error type:',error_type)
                time.sleep(1)
            print('\n    >>>> Error files saved in :'+str(errors_path)+'>>>>')      
        
        
        #Second optimization
        print('\nRefining previous structures at the :',levels[1].upper(),' level.\n')
        list_input_files_level2=os.listdir(path_level2)
        list_files=[]
        if len(list_input_files_level2)>=1:
            for file in list_input_files_level2:
                if file[-3:]=='com':
                    filename=os.path.join(path_level2,file)
                    list_files.append(filename)
            pool2=Pool(processes=4)
            result=pool2.map(RunGaussian,list_files)
            pool2.close()
            pool2.join()        
            
        #Identifying output files that converged 
        normal2_conv=[]
        error2_conv=[]
        for file in os.listdir(path_level2):
            if file[-3:]=='log':
                list_output_files_level2.append(file)
                filename=os.path.join(path_level2,file)
                INSTANCE=Opt_Freq(filename)                
                if INSTANCE.check_convergence()==1:
                    normal2_conv.append(filename)
                else:
                    error2_conv.append(filename)
                    shutil.copy(filename,errors_path2)
                    shutil.copy(filename[:-3]+'com',errors_path2)
                    os.remove(filename)
                    os.remove(filename[:-3]+'com')
        
        print('>>>> Reactants refined :',str(len(list_output_files_level2)))       
        
        if len(normal2_conv)>=1:
            print('>>>> Normally converged :',len(normal2_conv))   
            for filename in normal2_conv:
                print('    >>>>',filename, '\tSize(kB) :',int((os.stat(filename).st_size)/1024))
                time.sleep(1)
        
        if len(error2_conv)>=1:
            print('>>>> Error termination or job killed :',len(error2_conv))
            for filename in error2_conv:
                new_filename=os.path.join(errors_path2,os.path.split(filename)[1])
                error_type=get_error_type(new_filename)
                print('    >>>>',new_filename, '\tSize(kB) :',int((os.stat(new_filename).st_size)/1024),' Error type:',error_type)
                time.sleep(1)
            print('\n    >>>> Error files saved in :'+str(errors_path2)+'>>>>')      
            
            
        
    except BaseException as e:
        print('Error from the manager1 module (Optimize_R method)!!!')
        print(e)
    
    
    
    
    
def Optimize_C(my_path_1,levels,cwd):
    try:
        path_level2=os.path.join(my_path_1,'REFINED')
        errors_path=os.path.join(my_path_1,'ERROR_FILES')
        make_path(path_level2)
        make_path(errors_path)
        errors_path2=os.path.join(path_level2,'ERROR_FILES')
        make_path(errors_path2)
        list_input_files_level1=[]
        list_input_files_level2=[]
        list_output_files_level1=[]
        list_output_files_level2=[]
        
        if os.path.isdir(my_path_1):
            for file in os.listdir(my_path_1):
                if file[-3:]=='com':
                    filepath=os.path.join(my_path_1,file)
                    list_input_files_level1.append(filepath)
        
        print('\nGaussian called for calculation\nSystems: Cycloadduct(s)\nJob Type: OPT FREQ\nLEVEL: ',str(levels[0]).upper(),'\n')
        pool1=Pool(processes=4)
        result=pool1.map(RunGaussian,list_input_files_level1)
        pool1.close()
        pool1.join()   # The join() method block the execution until all the jobs fed into pool1 are terminated!!!      
        print('>>>> Gaussian Execution Terminated')
        
        #Identifying output files that converged and creating new inputs
        normal_conv=[]
        error_conv=[]
        for file in os.listdir(my_path_1):
            if file[-3:]=='log':
                list_output_files_level1.append(file)
                filename=os.path.join(my_path_1,file)
                INSTANCE=Opt_Freq(filename)                
                if INSTANCE.check_convergence()==1:
                    normal_conv.append(filename)
                    job_details='# opt freq '+str(levels[1])+' scf=(xqc, maxcycle=1024)'
                    INSTANCE.extract_gjf_file(path_level2,job_details)
                    shutil.copy(filename,path_level2)
                else:
                    error_conv.append(filename)
                    shutil.copy(filename,errors_path)
        
        print('>>>> Cycloadducts optimized :',str(len(list_output_files_level1)))       
        
        if len(normal_conv)>=1:
            print('>>>> Normally converged :',len(normal_conv))   
            for filename in normal_conv:
                print('    >>>>',filename, '\tSize(kB) :',int((os.stat(filename).st_size)/1024))
                time.sleep(1)
            print('\n    >>>> Structure(s) to be refined. Last geometries saved in :'+str(path_level2)+'>>>>')
        
        if len(error_conv)>=1:
            print('>>>> Error termination or job killed :',len(error_conv))
            for filename in error_conv:
                error_type=get_error_type(filename)
                print('    >>>>',filename, '\tSize(kB) :',int((os.stat(filename).st_size)/1024),' Error type:',error_type)
                time.sleep(1)
            print('\n    >>>> Error files saved in :'+str(errors_path)+'>>>>')      
        
        #Second optimization
        print('\nRefining previous structures at the :',levels[1].upper(),' level.\n')
        list_input_files_level2=os.listdir(path_level2)
        list_files=[]
        if len(list_input_files_level2)>=1:
            for file in list_input_files_level2:
                if file[-3:]=='com':
                    filename=os.path.join(path_level2,file)
                    list_files.append(filename)
            pool2=Pool(processes=4)
            result=pool2.map(RunGaussian,list_files)
            pool2.close()
            pool2.join()        
            
        #Identifying output files that converged 
        normal2_conv=[]
        error2_conv=[]
        for file in os.listdir(path_level2):
            if file[-3:]=='log':
                list_output_files_level2.append(file)
                filename=os.path.join(path_level2,file)
                INSTANCE=Opt_Freq(filename)                
                if INSTANCE.check_convergence()==1:
                    normal2_conv.append(filename)
                else:
                    error2_conv.append(filename)
                    shutil.copy(filename,errors_path2)
                    shutil.copy(filename[:-3]+'com',errors_path2)
                    os.remove(filename)
                    os.remove(filename[:-3]+'com')
        
        print('>>>> Cycloadduct(s) refined :',str(len(list_output_files_level2)))       
        
        if len(normal2_conv)>=1:
            print('>>>> Normally converged :',len(normal2_conv))   
            for filename in normal2_conv:
                print('    >>>>',filename, '\tSize(kB) :',int((os.stat(filename).st_size)/1024))
                time.sleep(1)
        
        if len(error2_conv)>=1:
            print('>>>> Error termination or job killed :',len(error2_conv))
            for filename in error2_conv:
                new_filename=os.path.join(errors_path2,os.path.split(filename)[1])
                error_type=get_error_type(new_filename)
                print('    >>>>',new_filename, '\tSize(kB) :',int((os.stat(new_filename).st_size)/1024),' Error type:',error_type)
                time.sleep(1)
            print('\n    >>>> Error files saved in :'+str(errors_path2)+'>>>>')      
            
            
            
        
    except BaseException as e:
        print('Error from the manager1 module (Optimize_C method)!!!')
        print(e)
    
    
    
def make_path(path):
    
    if os.path.isdir(path):
        shutil.rmtree(path)    
    os.mkdir(path)
        
def get_error_type(file_path):
    
    result=None
    
    ref_word1='Logic error in ASyTop'
    ref_word2='Atomic number out of range'
    ref_word3='Number of steps exceeded'
    ref_word4='Small interatomic distances encountered'
    ref_word5='Problem with the distance matrix'
    ref_word6='Atoms too close'
                   
    with open(file_path,'r') as myfile:            
        content=myfile.readlines()
        for line in content:
            if ref_word1 in line:
                result='Logic error in ASyTop'
                break
            if ref_word2 in line:
                result='Atomic number out of range'
                break                
            if ref_word3 in line:  
                result='Number of steps exceeded'
                break
            if ref_word4 in line or ref_word5 in line or ref_word6 in line:
                result='Problem with the distance matrix'
                break
            else:
                result='Unknown error'
                    

    return result             

    