#!/usr/bin/env python
# coding: utf-8

# In[1]:


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



from rdkit import Chem
import Geom_3D
import os
import sys
import time
import shutil
from Welcome import message
from configparser import ConfigParser
from manager1 import Optimize_R,Optimize_C
from manager2 import Optimize_TS
from manager3 import RunIRC
from pathconstr import Constr_irc_paths

message()

#####################################################################################
#                    Creating global variables  (job details)                       #
#####################################################################################

try:
    global cwd,D,levels_RC, levels_TS, cyclo, over_write,IRC_SS,IRC_LEVEL, SCRATCH,IRC_GEOMS_CONSTR
    cyclo=Chem.MolFromSmiles('C1=CCCCC1')
    cwd=os.getcwd()
    config=ConfigParser()
    config.read(os.path.join(cwd,"da.ini"))
    D=eval(config['job_details']['D_SPLITTING'])
    levels_RC=[x for x in (config['job_details']['CALC_LEVELS_RC']).split(';')]                      
    levels_TS=[x for x in (config['job_details']['CALC_LEVELS_TS']).split(';')]
    IRC_LEVEL=config['job_details']['CALC_LEVEL_IRC']
    RC_FLAG=int(config['flags']['RC_FLAG'])
    TS_FLAG=int(config['flags']['TS_FLAG'])
    IRC_FLAG=int(config['flags']['IRC_FLAG'])
    OVERWRITE_FLAG=int(config['flags']['OVERWRITE_FLAG'])
    IRC_GEOMS_CONSTR=int(config['flags']['IRC_GEOMS_CONSTR'])
    IRC_jobs=[int(x) for x in (config['job_details']['TS_ID_numbers']).split(',')]    
    NBR_IRC_POINT_PP=int(config['job_details']['NBR_IRC_POINT_PP'])  
    IRC_SS = int(config['job_details']['IRC_STEP_SIZE'])  
    NBR_PATHS=int(config['job_details']['NBR_PATHS'])
    SCRATCH =int(config['flags']['SCRATCH'])
    

except BaseException as e:    
    sys.exit('Possible error in the da.ini file.\nThe execution of AMADAR has been aborted')


####################################################################################
#      checking the existence of the SMILES.txt file  and making a copy of it      #
####################################################################################
test_dir=os.path.join(cwd,'test')
if os.path.isdir(test_dir):
    shutil.rmtree(test_dir)
os.mkdir(test_dir)

if os.path.isfile(os.path.join(cwd,'SMILES.txt')):
    shutil.copy(os.path.join(cwd,'SMILES.txt'),os.path.join(cwd,'SMILES_temp.txt'))
    shutil.copy(os.path.join(cwd,'SMILES.txt'),test_dir)
else:
    sys.exit("SMILES.txt does not exist.\nFor this calculation we require a SMILES.txt with a series of smiles strings for DA adducts")

####################################################################################
#               creating the errors.txt and the AOI.csv files                      #
####################################################################################

    
if os.path.isfile(os.path.join(cwd,'errors.txt')):
    os.remove(os.path.join(cwd,'errors.txt'))
else:
    error_file=open('errors_noted.txt','a')
    error_file.close()
    
Geom_3D.Gen_AOI_file(cwd)
    
####################################################################################
#                          INITIALIZATIONN                                         #
#              Generation of input 3D geometries in tree steps:                    #
#                  1. Creation of RDKit mol objects                                #
#                  2. Retro-DielsAlder and identification of the reactive site     #
#                  3. Emdedding, conformational search and UFF optimization        #
####################################################################################



def main():
    
    try: 
        
        smiles_path=os.path.join(cwd,'SMILES_temp.txt')
        if os.path.isfile(smiles_path):        
            with open(smiles_path,'r') as file:        
                content=file.readlines()            
            if len(content)!=0 and content[0]!='\n':            
                print('SMILES strings uploaded!!!\n') 
                
                ####################################################################################
                #      Checking the SMILES strings to make sure they are real DA cycloadducts      #
                ####################################################################################
                cycloadducts_pos=[]
                not_cycloaddcuts=[]
                for i in range(len(content)):
                    if Chem.MolFromSmiles(content[i]).HasSubstructMatch(cyclo)==True:
                        cycloadducts_pos.append(i)
                    else:
                        not_cycloaddcuts_pos.append(i)
                
                if len(cycloadducts_pos)>=1:
                    time.sleep(2)
                    print('    <<<<',len(cycloadducts_pos),' DA cycloadducts identified in the SMILES.txt file')             
                    path_reagents=os.path.join(cwd,'R')
                    path_cycloadd=os.path.join(cwd,'C')
                    path_pgts=os.path.join(cwd,'TS')    
                    paths=[path_reagents,path_cycloadd, path_pgts]    
                    for path in paths:
                        if os.path.isdir(path):
                            if OVERWRITE_FLAG==1:
                                shutil.rmtree(path)                        
                                os.mkdir(path)
                        else:
                            os.mkdir(path)                       
                                                      
                    time.sleep(2)
                    if SCRATCH==1:
                        sms2='\nGeneration of initial 3D geometries in progress\n'    
                        print(sms2)
                        for i in cycloadducts_pos:       
                            mol=Chem.MolFromSmiles(content[i])               
                            Geom_3D.Gen_gjf_file_reagents(mol,path_reagents,i+1,levels_RC[0])  
                            Geom_3D.Gen_gjf_file_adducts(mol,path_cycloadd,i+1, levels_RC[0])            
                            Geom_3D.Gen_gjf_file_ts(mol,path_pgts,i+1,levels_TS[0],D)            
                            if (i+1)%100==0:                
                                print('%s reactions treated' %(i+1))     
                        time.sleep(1)
                        print('    <<<< Input files successfully generated with 3D geometries obtained at UFF level >>>')
                        time.sleep(2) 
                    else:
                        print('\n    <<<< Input generation step skipped!!!')
                        time.sleep(2)
                        print('\n    <<<< Make sure at least one of the the R, C or TS directories already contain files to process')
                        time.sleep(2)

                    ####################################################################################
                    #       Checking job flags before running electronic structure calculations        #
                    ####################################################################################


                    if RC_FLAG==1:
                        Optimize_R (path_reagents,levels_RC,cwd)
                        Optimize_C(path_cycloadd,levels_RC,cwd)
                    if TS_FLAG==1 :
                        Optimize_TS(path_pgts,levels_TS,cwd)

                    if IRC_FLAG==1 : 

                        if isinstance(IRC_jobs,list) and isinstance(NBR_IRC_POINT_PP,int):

                            RunIRC(cwd,NBR_IRC_POINT_PP,IRC_jobs,IRC_SS,IRC_LEVEL)
                        else:
                            print('IRC_job aborted!!!')
                            print('Make sur IRC_jobs = list and NBR_IRC_POINT_PP=int in the da.ini file')
                            
                    if IRC_GEOMS_CONSTR==1:
                        Constr_irc_paths(cwd)
                
                else:
                    
                    sys.exit('Not Diels-Alder cycloadduct found in the SMILES.txt file\nJob has been aborted')                    
        
            else:
                print('SMILES.txt file uploaded, but may be empty or has a blank line at the top. Please check it.\nProcess killed!!!')
                
    except BaseException as e:            
        error_message=str(e)
        sys.exit(error_message)         
        
    #Removing the temporary file
    if os.path.isfile(os.path.join(cwd,'SMILES_temp.txt')):
        os.remove(os.path.join(cwd,'SMILES_temp.txt'))  
    #Returning to the root directory
    os.chdir(cwd)
    
if __name__=='__main__':

    main()


# In[ ]:




