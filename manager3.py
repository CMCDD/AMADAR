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
from TS import Opt_TS
from IRC import create_IRC_input
from Gaussian import RunGaussian
from multiprocessing import Pool
from Min import Opt_Freq
import shutil
import time
from configparser import ConfigParser




def RunIRC(cwd,Nbr_IRC_pp,IRC_jobs,step_size,IRC_LEVEL):
    
    config=ConfigParser()
    config.read(os.path.join(cwd,"da.ini"))
    OVERWRITE_FLAG=int(config['flags']['OVERWRITE_FLAG'])
    
    print('\nChecking the predicted TS(s)')
    list_ts=[]
    list_ts_located=[]
    list_ts_not_located=[]
    list_output_files_irc=[]
    try:
        ts_path=os.path.join(cwd,'TS/GTS/TS')
        irc_path=os.path.join(cwd,'IRC')
        make_dir(irc_path,OVERWRITE_FLAG)
        ts_not_located=os.path.join(irc_path,'TS_NOT_FOUND')
        make_dir(ts_not_located,OVERWRITE_FLAG)
        irc_errors=os.path.join(irc_path,'ERROR_FILES')
        make_dir(irc_errors,OVERWRITE_FLAG)
        
        if os.path.isdir(ts_path):            
            for file in os.listdir(ts_path):
                if IRC_jobs[0].upper()=='ALL':
                    if file[-3:]=='log':
                        filename=os.path.join(ts_path,file)
                        list_ts.append(filename)
                else:
                    if file[-3:]=='log' and int(file[4:-6]) in IRC_jobs:
                        print(file)
                        filename=os.path.join(ts_path,file)
                        list_ts.append(filename)
                    
            #Checking and reporting TS location 
            for filename in list_ts:            
                INSTANCE=Opt_TS(filename)
                Vibr_anal=INSTANCE.Analyze_vibr_freq()
                imag_freq=Vibr_anal[0]
                pseudo_imag_freq=Vibr_anal[1]
                if (len(imag_freq)+len(pseudo_imag_freq))==1:
                    list_ts_located.append(filename)
                else:
                    list_ts_not_located.append(filename)
                    shutil.copy(filename,ts_not_located)
                    shutil.copy(filename[:-3]+'com',ts_not_located)
                    os.remove(filename)
                    os.remove(filename[:-3]+'com')

            if len(list_ts_located)>=1:
                print('>>>> TS(s) located :',len(list_ts_located))   
                for filename in list_ts_located:
                    imag_freq=Opt_TS(filename).Analyze_vibr_freq()[0]
                    print('    >>>>',filename, '\tSize(kB) :',int((os.stat(filename).st_size)/1024),'\tImagFreq :',imag_freq[0])
                    time.sleep(1)

            if len(list_ts_not_located)>=1:
                print('>>>> Error termination or job killed :',len(list_ts_not_located))
                for filename in list_ts_not_located:
                    new_filename=os.path.join(ts_not_located,os.path.split(filename)[1])
                    error_type=get_error_type(new_filename)
                    print('    >>>>',new_filename, '\tSize(kB) :',int((os.stat(new_filename).st_size)/1024),' Error type:',error_type)
                    time.sleep(1)
                print('\n    >>>> TS(s) not located. Files saved in :'+str(ts_not_located)+'>>>>')        
            time.sleep(2)


            #Creating inputs for IRC calculations
            for filename in list_ts_located:
                INSTANCE2=create_IRC_input(filename)
                job_details='# irc=(maxpoints='+str(Nbr_IRC_pp)+',recorrect=test,calcall,stepsize='+str(step_size)+') '+str(IRC_LEVEL)
                INSTANCE2.get_unique_input(irc_path,job_details)

            time.sleep(2)
            print('\n    >>>> '+str(len(list_ts_located))+' IRC calculation input(s) created')

            #Running IRC calculations
            list_files=[]
            for file in os.listdir(irc_path):
                if file[-3:]=='gjf' or file[-3:]=='com':
                    filename=os.path.join(irc_path,file)
                    list_files.append(filename)   #List of input files ofr IRC calculations
            print('\nRunning IRC calculations')
            print('This step may take so long. Grab your cup of tea and wait patiently.\nTalk to you soon.')
            pool=Pool(processes=4)
            result=pool.map(RunGaussian,list_files)
            pool.close()
            pool.join()
            print('\n>>>> Gaussian execution terminated')
            
            #Checking IRC calculations that converged
            list_irc_norm_conv=[]
            list_irc_error_conv=[]
            for file in os.listdir(irc_path):
                if file[-3:].upper()=='LOG':
                    filename=os.path.join(irc_path,file) 
                    list_output_files_irc.append(filename)                               
                    if Opt_Freq(filename).check_convergence()==1:
                        print(filename)
                        list_irc_norm_conv.append(filename)
                    else:
                        list_irc_error_conv.append(filename)
                        shutil.copy(filename,irc_errors)
                        shutil.copy(filename[:-3]+'gjf',irc_errors)
                        os.remove(filename)
                        os.remove(filename[:-3]+'gjf')

            print('>>>> IRC calculations completed :',str(len(list_output_files_irc)))

            if len(list_irc_norm_conv)>=1:
                print('>>>> Normally converged :',len(list_irc_norm_conv))   
                for filename in list_irc_norm_conv:
                    print('    >>>>',filename, '\tSize(kB) :',int((os.stat(filename).st_size)/1024))
                    time.sleep(1)

            if len(list_irc_error_conv)>=1:
                print('>>>> Error termination or job killed :',len(list_irc_error_conv))
                for filename in list_irc_error_conv:
                    new_filename=os.path.join(irc_errors,os.path.split(filename)[1])
                    error_type=get_error_type(new_filename)
                    print('    >>>>',new_filename, '\tSize(kB) :',int((os.stat(new_filename).st_size)/1024),' Error type:',error_type)
                    time.sleep(1)
                print('\n    >>>> Error files saved in :'+str(irc_errors)+'>>>>')             
            
        
        
        else:
            print('TS path not found')    
        
        
    
    except BaseException as e:
        
        print('Error in myIRC module')
        print(e)


def make_dir(path,OVERWRITE_FLAG):    
    
    if os.path.isdir(path)==False:
        os.mkdir(path)
    else:
        if OVERWRITE_FLAG==1:
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

    
