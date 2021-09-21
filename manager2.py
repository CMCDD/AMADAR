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
from Guess import RCO_guess_ts
from TS import Opt_TS
import shutil
import time




def Optimize_TS(my_path_1,levels,root):
    
    try:
        global path_level2,path_level3,paths
        path_gts=os.path.join(my_path_1,'GTS')
        errors_pgts=os.path.join(my_path_1,'ERROR_FILES')
        for path in [path_gts,errors_pgts]:
            make_path(path)    
        path_ts=os.path.join(path_gts,'TS')
        errors_gts=os.path.join(path_gts,'ERROR_FILES')
        for path in [path_ts,errors_gts]:
            make_path(path)            
        errors_ts=os.path.join(path_ts,'ERROR_FILES')
        make_path(errors_ts)
        
        list_input_files_pgts=[]
        list_input_files_gts=[]
        list_input_files_ts=[]
        list_output_files_pgts=[]
        list_output_files_gts=[]
        list_output_files_ts=[]
        
        #Identifying input files for pgts calcs
        if os.path.isdir(my_path_1):
            for file in os.listdir(my_path_1):
                if file[-3:]=='com':
                    filepath=os.path.join(my_path_1,file)
                    list_input_files_pgts.append(filepath)
        else:
            print('Path not found!!!')
            
        print('\nGaussian called for calculation\nSystems: Pseudo-guess TS(s)\nJob Type: RCO\nLEVEL: ',str(levels[0]).upper(),'\n')
        pool1=Pool(processes=4)
        result=pool1.map(RunGaussian,list_input_files_pgts)
        pool1.close()
        pool1.join()   # The join() method block the execution until all the jobs fed into pool1 are terminated!!!      
        print('>>>> Gaussian Execution Terminated')
        
        #Identifying output files that converged and creating new inputs
        normal_conv_pgts=[]
        error_conv_pgts=[]
        for file in os.listdir(my_path_1):
            if file[-3:]=='log':
                list_output_files_pgts.append(file)
                filename=os.path.join(my_path_1,file)
                INSTANCE=RCO_guess_ts(filename)                
                opt_steps=INSTANCE.get_opt_steps_energies() 
                list_conf_energies=opt_steps[0]
                target_conf_energy=opt_steps[1]
                target_conf_index=(list_conf_energies.index(target_conf_energy)+1)            
                if Opt_Freq(filename).check_convergence()==1:
                    normal_conv_pgts.append(filename)
                    job_details='# opt=(calcall,ts,noeigentest,maxstep=5) '+str(levels[1])+' scf=xqc'                    
                    INSTANCE.extract_last_conf(path_gts,target_conf_index,job_details)
                    #print(filename)
                else:
                    error_conv_pgts.append(filename)
                    shutil.copy(filename,errors_pgts)
                    shutil.copy(filename[:-3]+'com',errors_pgts)
                    os.remove(filename)
                    os.remove(filename[:-3])
        
        
        print('>>>> Pseudo-guess TSs optimized :',str(len(list_output_files_pgts)),'\t',time.ctime())
        
        #Checking optimization results
        if len(normal_conv_pgts)>=1:
            print('>>>> Normally converged :',len(normal_conv_pgts))   
            for filename in normal_conv_pgts:
                print('    >>>>',filename, '\tSize(kB) :',int((os.stat(filename).st_size)/1024))
                time.sleep(1)
            print('\n    >>>> Highest energy configuration(s) saved in :'+str(path_gts)+'>>>>')
        
        if len(error_conv_pgts)>=1:
            print('>>>> Error termination or job killed :',len(error_conv_pgts))
            for filename in error_conv_pgts:
                new_filename=os.path.join(errors_pgts,os.path.split(filename)[1])
                error_type=get_error_type(new_filename)
                print('    >>>>',new_filename, '\tSize(kB) :',int((os.stat(new_filename).st_size)/1024),' Error type:',error_type)
                time.sleep(1)
            print('\n    >>>> Error files saved in :'+str(errors_pgts)+'>>>>')              
        
        #Relaxing pguess towards the guess TS
        print('\nRelaxing pseudo-guess TS(s) at the ',levels[1].upper(),' level.\n')
        
        list_input_files_gts=os.listdir(path_gts)
        list_files=[]
        if len(list_input_files_gts)>=1:
            for file in list_input_files_gts:
                if file[-3:]=='com':
                    filename=os.path.join(path_gts,file)
                    list_files.append(filename)
            pool2=Pool(processes=4)
            result=pool2.map(RunGaussian,list_files)
            pool2.close()
            pool2.join() 
            
        #Identifying guess TSs that converged and constructing new inputs for the TS 
        normal_conv_gts=[]
        error_conv_gts=[]
        for file in os.listdir(path_gts):
            if file[-3:]=='log':
                list_output_files_gts.append(file)
                filename=os.path.join(path_gts,file)            
                if Opt_Freq(filename).check_convergence()==1:
                    normal_conv_gts.append(filename)
                    
                    
                else:
                    error_conv_gts.append(filename)
                    shutil.copy(filename,errors_gts)
                    shutil.copy(filename[:-3]+'com',errors_gts)
                    os.remove(filename)
                    os.remove(filename[:-3]+'com')
                    
        print('>>>> Guess TSs predicted :',str(len(list_output_files_gts)), '\t', time.ctime())       
        
        if len(normal_conv_gts)>=1:
            print('>>>> Normally converged :',len(normal_conv_gts))   
            for filename in normal_conv_gts:
                print('    >>>>',filename, '\tSize(kB) :',int((os.stat(filename).st_size)/1024))
                time.sleep(1)
        
        if len(error_conv_gts)>=1:
            print('>>>> Error termination or job killed :',len(error_conv_gts))
            for filename in error_conv_gts:
                new_filename=os.path.join(errors_gts,os.path.split(filename)[1])
                error_type=get_error_type(new_filename)
                print('    >>>>',new_filename, '\tSize(kB) :',int((os.stat(new_filename).st_size)/1024),' Error type:',error_type)
                time.sleep(1)
            print('\n    >>>> Error files saved in :'+str(errors_path_gts)+'>>>>')
            
        #New relaxtion of GTS into the likely TS. 
        print('\nNew relaxtion of guess TSs into the likely TS at the ',levels[2].upper(),' level.\n')
        for file in os.listdir(path_gts):
            if file[-3:]=='log':
                filename=os.path.join(path_gts,file)
                INSTANCE=Opt_TS(filename)
                job_details='# opt=(calcall,ts,noeigentest) freq '+str(levels[2])+' scf=xqc'
                INSTANCE.extract_ts_gjf_file(path_ts,job_details)
        
        list_input_files_ts=os.listdir(path_ts)
        list_files=[]
        if len(list_input_files_ts)>=1:
            for file in list_input_files_ts:
                if file[-3:]=='com':          #Not necessary but in case someone volontary insert a non .com file in the folder
                    filename=os.path.join(path_ts,file)
                    list_files.append(filename)
            pool2=Pool(processes=4)
            result=pool2.map(RunGaussian,list_files)
            pool2.close()
            pool2.join() 
        
        #Identifying TSs calcs that converged 
        normal_conv_ts=[]
        error_conv_ts=[]
        for file in os.listdir(path_ts):
            if file[-3:]=='log':
                list_output_files_ts.append(file)
                filename=os.path.join(path_ts,file)            
                if Opt_Freq(filename).check_convergence()==1:
                    normal_conv_ts.append(filename)
                else:
                    error_conv_ts.append(filename)
                    shutil.copy(filename,errors_ts)
                    shutil.copy(filename[:-3]+'com',errors_ts)
                    os.remove(filename)
                    os.remove(filename[:-3]+'com')
        
        print('>>>> TSs predicted :',str(len(list_output_files_ts)))
        
        if len(normal_conv_ts)>=1:
            print('>>>> Normally converged :',len(normal_conv_ts))   
            for filename in normal_conv_ts:
                print('    >>>>',filename, '\tSize(kB) :',int((os.stat(filename).st_size)/1024))
                time.sleep(1)
        
        if len(error_conv_ts)>=1:
            print('>>>> Error termination or job killed :',len(error_conv_ts))
            for filename in error_conv_ts:
                new_filename=os.path.join(errors_ts,os.path.split(filename)[1])
                error_type=get_error_type(new_filename)
                print('    >>>>',new_filename, '\tSize(kB) :',int((os.stat(new_filename).st_size)/1024),' Error type:',error_type)
                time.sleep(1)
            print('\n    >>>> Error files saved in :'+str(errors_ts)+'>>>>') 
            
        #Last step checking stationary points
        #Checking_TS(path_ts)
        #Correct error files 
       
    except BaseException as e:
        
        print('Error in the manager2 module')
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

    
        