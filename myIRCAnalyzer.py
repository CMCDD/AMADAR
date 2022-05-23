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



import os
import shutil
from Welcome import message
from RFA import RFA_1file
from RFA import MultiReactionForceConstant1F,MultiReactionForce1F,MultiPotentialEnergy1F
from RFD import RF_atomic_decomp,RFC_atomic_decomp
from WBOA import WBO_analysis_1F
import time
from configparser import ConfigParser
import sys


cwd=os.getcwd()
message()
print('Welcome to the IRC path analysis module\n')

def main():   
    
    global irc_folder
    irc_folder=os.path.join(cwd,'IRC')
    
    config=ConfigParser()
    config.read(os.path.join(cwd,"da.ini"))
    RFA_FLAG=int(config['flags']['RFA_FLAG'])
    RFD_FLAG=int(config['flags']['RFD_FLAG'])
    WBOA_FLAG=int(config['flags']['WBOA_FLAG'])
    
    config2=ConfigParser()
    config2.read(os.path.join(cwd,"analysis.ini"))
    RFA_IDs = [int(x) for x in (config2['RFA']['Unq_RFA']).split(',')]
    Multiple_RF = [int(x) for x in (config2['RFA']['Multiple_RF']).split(',')]
    Multiple_RE = [int(x) for x in (config2['RFA']['Multiple_RE']).split(',')]
    Multiple_RFC =[int(x) for x in (config2['RFA']['Multiple_RFC']).split(',')]    
    RFD_JOB_ID = int(config2['RFD']['JOB_ID'])
    RFD_ATOMS = [int(x) for x in (config2['RFD']['ATOMS']).split(',')]
    RFD_FRAG = [x for x in (config2['RFD']['FRAG']).split(';')]
    RFD_FRAG_NAMES = [x for x in (config2['RFD']['FRAG_NAMES']).split(',')]
    WBOA_ID = [int(x) for x in (config2['IRC_PATHS']['WBOA_ID']).split(',')]
    
    if RFA_FLAG==1:        
        
        time.sleep(1)
        RFA_IDs = [int(x) for x in (config2['RFA']['Unq_RFA']).split(',')]
        Multiple_RF = [int(x) for x in (config2['RFA']['Multiple_RF']).split(',')]
        Multiple_RE = [int(x) for x in (config2['RFA']['Multiple_RE']).split(',')]
        Multiple_RFC =[int(x) for x in (config2['RFA']['Multiple_RFC']).split(',')] 
        
        ###################################################
        if len(RFA_IDs)==1 and RFA_IDs[0]==-1:
            print('\nReaction force analysis will be run for all the files in the IRC directory')
            time.sleep(2)
            list_paths=[]
            for file in os.listdir(irc_folder):
                if file.endswith('.log'):
                    filename=os.path.join(irc_folder,file)
                    list_paths.append(filename)
            count=1        
            for path in list_paths:
                if os.path.isfile(path):
                    if check_IRC_validity(path)=='RELIABLE':
                        path_rfa=RunRFA(path,count)     
                    else:
                        print(path,' : IRC path not reliable')
                        print('Possible solutions : ')
                        print('1. Regenerate TS from the very begining with a different D value before rerunning IRC calculation') 
                        print('2. Determine the IRC path using the 2files approach')
                        print('For the latter option, set NBR_PATHS=2 in the da.ini file')
                else:
                    print(path,' not found')
                count+=1   
                
            if os.path.isfile(os.path.join(path_rfa,'rfa.txt')):
                file=open(os.path.join(path_rfa,'rfa.txt'),'a')
                file.write('\nJob terminated successfully\n')
                file.write(time.ctime())                
        
            
        elif len(RFA_IDs)>=1 and RFA_IDs[0]!=0:
            list_paths=[]
            for file in os.listdir(irc_folder):
                if file.endswith('.log') and int(file[5:-6]) in RFA_IDs:
                    filename=os.path.join(irc_folder,file)
                    list_paths.append(filename)
            count=1        
            for path in list_paths:
                if os.path.isfile(path):
                    if check_IRC_validity(path)=='RELIABLE':
                        path_rfa=RunRFA(path,count)     
                    else:
                        print(path,' : IRC path not reliable')
                        print('Possible solutions : ')
                        print('1. Regenerate TS from the very begining with a different D value before rerunning IRC calculation') 
                        print('2. Determine the IRC path using the 2files approach')
                        print('For the latter option, set NBR_PATHS=2 in the da.ini file')
                else:
                    print(path,' not found')
                count+=1   
                
            if os.path.isfile(os.path.join(path_rfa,'rfa.txt')):
                file=open(os.path.join(path_rfa,'rfa.txt'),'a')
                file.write('\nJob terminated successfully\n')
                file.write(time.ctime())
        else:
            sys.exit('Please make sure the IRC_IDs provided correspond to any file in the IRC directory!!!')
                
        ###################################################
                
        if len(Multiple_RF)>=1 and Multiple_RF[0]!=0:
            list_paths=[]
            for file in os.listdir(irc_folder):
                if file.endswith('.log') and int(file[5:-6]) in Multiple_RF:
                    filename=os.path.join(irc_folder,file)
                    list_paths.append(filename)
            
            new_list_paths=[]        
            for path in list_paths:
                if os.path.isfile(path):
                    if check_IRC_validity(path)=='RELIABLE':
                        new_list_paths.append(path)
                    else:
                        print(path,' : IRC path not reliable')
                        print('Possible solutions : ')
                        print('1. Regenerate TS from the very begining with a different D value before rerunning IRC calculation') 
                        print('2. Determine the IRC path using the 2files approach')
                        print('For the latter option, set NBR_PATHS=2 in the da.ini file')
                
                else:
                    print(path,' not found')
            RunMultiple_RF(new_list_paths)
            
        ###################################################
        
        if len(Multiple_RE)>=1 and Multiple_RE[0]!=0:
            list_paths=[]
            for file in os.listdir(irc_folder):
                if file.endswith('.log') and int(file[5:-6]) in Multiple_RE:
                    filename=os.path.join(irc_folder,file)
                    list_paths.append(filename)
                    
            new_list_paths=[]    
            for path in list_paths:
                if os.path.isfile(path):
                    if check_IRC_validity(path)=='RELIABLE':
                        new_list_paths.append(path)
                    else:
                        print(path,' : IRC path not reliable')
                        print('Possible solutions : ')
                        print('1. Regenerate TS from the very begining with a different D value before rerunning IRC calculation') 
                        print('2. Determine the IRC path using the 2files approach')
                        print('For the latter option, set NBR_PATHS=2 in the da.ini file')
                else:
                    print(path,' not found')
                    
            RunMultiple_RE(new_list_paths)
            
        ###################################################
            
        if len(Multiple_RFC)>=1 and Multiple_RFC[0]!=0:
            list_paths=[]
            for file in os.listdir(irc_folder):
                if file.endswith('.log') and int(file[5:-6]) in Multiple_RFC:
                    filename=os.path.join(irc_folder,file)
                    list_paths.append(filename)
                    
            new_list_paths=[]
            for path in list_paths:
                if os.path.isfile(path):
                    if check_IRC_validity(path)=='RELIABLE':
                        new_list_paths.append(path)
                    else:
                        print(path,' : IRC path not reliable')
                        print('Possible solutions : ')
                        print('1. Regenerate TS from the very begining with a different D value before rerunning IRC calculation') 
                        print('2. Determine the IRC path using the 2files approach')
                        print('For the latter option, set NBR_PATHS=2 in the da.ini file')
                else:
                    print(path,' not found')
            RunMultiple_RFC(new_list_paths)
            
        ###################################################
            
        ###################################################
            
        ###################################################
            
    if RFD_FLAG==1: 
            
        for file in os.listdir(irc_folder):
            if file.endswith('.log') and int(file[5:-6])==RFD_JOB_ID:
                if len(RFD_ATOMS)>=1 and RFD_ATOMS[0]!=0:
                    if check_IRC_validity(os.path.join(irc_folder,file))=='RELIABLE':
                        filename=os.path.join(irc_folder,file)
                        Run_RF_AtomicDecomposition(filename,RFD_ATOMS)
                        Run_RFC_AtomicDecomposition(filename,RFD_ATOMS)
                    else:
                        print(os.path.join(irc_folder,file),' : IRC path not reliable')
                        print('Possible solutions : ')
                        print('1. Regenerate TS from the very begining with a different D value before rerunning IRC calculation') 
                        print('2. Determine the IRC path using the 2files approach')
                        print('For the latter option, set NBR_PATHS=2 in the da.ini file')
            
        ###################################################################################################
        
            
        ###################################################################################################    
        if len(RFD_FRAG)>=1 and RFD_FRAG[0]!=0:  
            for file in os.listdir(irc_folder):
                if file.endswith('.log') and int(file[5:-6])==RFD_JOB_ID:
                    if check_IRC_validity(os.path.join(irc_folder,file))=='RELIABLE':
                        if len(RFD_FRAG)==len(RFD_FRAG_NAMES):
                            filename=os.path.join(irc_folder,file)
                            list_fragments=[]
                            list_fragment_names=[]
                            for i in range(len(RFD_FRAG)):
                                list_fragments.append([])
                                list_fragment_names.append(RFD_FRAG_NAMES[i])
                                for j in RFD_FRAG[i].split(','):
                                    index=int(j)
                                    list_fragments[i].append(index)
                    
                            Run_RF_Fragment_Decomposition(filename,list_fragments,list_fragment_names)
                            Run_RFC_Fragment_Decomposition(filename,list_fragments,list_fragment_names)                   
                        
                        
                        else:
                            print('Number of frags != Frag names ')
                            print('Check the analysis.ini file and make sure there are as many frags as their names')
                            
                    else:
                        print(os.path.join(irc_folder,file),' : IRC path not reliable')
                        print('Possible solutions : ')
                        print('1. Regenerate TS from the very begining with a different D value before rerunning IRC calculation') 
                        print('2. Determine the IRC path using the 2files approach')
                        print('For the latter option, set NBR_PATHS=2 in the da.ini file')           
                    
                                       
        ###################################################################################################  
    if WBOA_FLAG==1:
        
        if len(WBOA_ID)>=1 and WBOA_ID[0]!=0:
        
            for ID in WBOA_ID:
            
                Run_WBOAnalysis(cwd,ID)
            
            
        ###################################################################################################  
    
def RunRFA(path,count):
    
    try:
        
        path_rfa=os.path.join(cwd,'rfa')
        make_dir2(path_rfa)
        print('\nReaction Force Analysis (RFA) in progress...\n')
        time.sleep(2)
        filename=path
        
        if count==1:
            if os.path.isfile(os.path.join(path_rfa,'rfa.txt')):
                os.remove(os.path.join(path_rfa,'rfa.txt'))
                myfile=open(os.path.join(path_rfa,'rfa.txt'),'a')
                myfile.write('Welcome to the AMADAR output file\n\n')
                myfile.write('COMPUTATIONAL MECHANISTIC CHEMISTRY AND DRUG DISCOVERY\n')
                myfile.write('Rhodes University\n')
                myfile.write('Department of Chemistry\n\n')
                myfile.write('TASK PERFORMED : REACTION FORCE ANALYSIS\n')
                myfile.write(time.ctime()+'\n')
                myfile.write('+-'*45)
                myfile.write('\n')
            else:
                myfile=open(os.path.join(path_rfa,'rfa.txt'),'a')
                myfile.write('Welcome to the AMADAR output file\n\n')
                myfile.write('COMPUTATIONAL MECHANISTIC CHEMISTRY AND DRUG DISCOVERY\n')
                myfile.write('Rhodes University\n')
                myfile.write('Department of Chemistry\n\n')
                myfile.write('TASK PERFORMED : REACTION FORCE ANALYSIS\n')
                myfile.write(time.ctime()+'\n')
                myfile.write('+-'*45)
                myfile.write('\n')
                
        else:
            myfile=open(os.path.join(path_rfa,'rfa.txt'),'a')
            
        if os.path.split(filename)[1][-3:].upper()=='LOG':
            instance=RFA_1file(filename)
            RFA=instance.plot(path_rfa)
            En_dec=instance.dec_activ_energy()
            EA,RE,TSE=RFA[0],RFA[1],RFA[2]
            EA_ALPHA_R=En_dec[0]
            EA_TS_ALPHA=En_dec[1]
            EA_GAMMA_TS=En_dec[2] 
            EA_P_GAMMA=En_dec[3]  
            RFC_MIN_POS=En_dec[4]
            print(En_dec)
            file=os.path.split(filename)[1]   
            print('>>>',file[3:-4])
            print('    >>>> Activation energy (kcal/mol) =',round(EA,2))
            time.sleep(0.40)
            print('    >>>> Reaction energy (kcal/mol) =',round(RE,2))
            time.sleep(0.40)
            print('    >>>> TS extent (ξ) =',round(TSE,2))
            time.sleep(0.40)
            print('    >>>> React.Forc.Const_MIN_po =',round(RFC_MIN_POS,2))
            time.sleep(0.40)
            print('    >>>> Pot. energy (R <-> ALPHA) (kcal/mol) =',round(EA_ALPHA_R,2))
            time.sleep(0.40)
            print('    >>>> Pot. energy (ALPHA <-> TS) (kcal/mol) =',round(EA_TS_ALPHA,2))
            time.sleep(0.40)
            print('    >>>> Pot. energy (GAMMA <-> TS) (kcal/mol) =',round(EA_GAMMA_TS,2))
            time.sleep(0.40)
            print('    >>>> Pot. energy (P <-> GAMMA) (kcal/mol) =',round(EA_P_GAMMA,2))
            
            time.sleep(0.40)
            myfile.write('\n\n SYSTEM : '+file[3:-4]+'\n')
            myfile.write('\n\n')
            myfile.writelines(['Activation energy (kcal/mol) =',' ',str(round(EA,2)),'\n'])
            myfile.writelines(['Reaction energy (kcal/mol) =',' ',str(round(RE,2)),'\n'])
            myfile.writelines(['TS extent (ξ) =','\t',str(round(TSE,2)),'\n'])
            myfile.writelines(['Pot. energy (R <-> ALPHA) (kcal/mol) =',' ',str(round(EA_ALPHA_R,2)),'\n'])
            myfile.writelines(['Pot. energy (ALPHA <-> TS) (kcal/mol) =',' ',str(round(EA_TS_ALPHA,2)),'\n'])
            myfile.writelines(['Pot. energy (GAMMA <-> TS) (kcal/mol) =',' ',str(round(EA_GAMMA_TS,2)),'\n'])
            myfile.writelines(['Pot. energy (P <-> GAMMA) (kcal/mol) =',' ',str(round(EA_P_GAMMA,2)),'\n'])  
            myfile.close()
            
    except BaseException as e:
        print(e)
        print('Error in the myIRCAnalyzer module')
        file=open('errors_noted.txt','a')
        file.write(str(e)+str(time.ctime()))
        file.close()
        
    return path_rfa

def RunMultiple_RFC(list_paths):
    
    try:
        
        path_rfa=os.path.join(cwd,'rfc_multi')
        make_dir2(path_rfa)
        print('\nMultipath analyzer launched. Reaction force constant curves to be plotted... \n')
        time.sleep(2)
        list_irc=[]
        for filename in list_paths:
            if os.path.split(filename)[1][-3:].upper()=='LOG':
                list_irc.append(filename)        
        instance=MultiReactionForceConstant1F(list_irc)
        instance.plot(os.path.join(path_rfa,'mrfc.png'),[7,5])   #[5,3] is the fig size
        
    except BaseException as e:
        
        print(e)

def RunMultiple_RF(list_paths):
    
    try:
        
        
        path_rfa=os.path.join(cwd,'rf_multi')
        make_dir2(path_rfa)
        print('\nMultipath analyzer launched. Reaction force curves to be plotted... \n')
        time.sleep(2)
        list_irc=[]
        for filename in list_paths:
            if os.path.split(filename)[1][-3:].upper()=='LOG':
                list_irc.append(filename)   
                
        instance=MultiReactionForce1F(list_irc)
        instance.plot(os.path.join(path_rfa,'mrf.png'),[7,5])   #[5,3] is the fig size
        
    except BaseException as e:
        
        print(e)

def RunMultiple_RE(list_paths):
    
    try:
        
        path_rfa=os.path.join(cwd,'re_multi')
        make_dir2(path_rfa)
        print('\nMultipath analyzer launched. Reaction energy curves to be plotted... \n')
        time.sleep(2)
        list_irc=[]
        for filename in list_paths:
            if os.path.split(filename)[1][-3:].upper()=='LOG':
                list_irc.append(filename)
        
        instance=MultiPotentialEnergy1F(list_irc)
        instance.plot(os.path.join(path_rfa,'mre.png'),[7,5])   #[5,3] is the fig size
        
    except BaseException as e:
        
        print(e)
    
def Run_RF_AtomicDecomposition(path,list_atomic_index):
    
    try:
        
        path_rfd=os.path.join(cwd,'rfd')
        make_dir2(path_rfd)
        
        instance=RF_atomic_decomp(path)  
        print('   <<<< Processing...')
        instance.plot_multiple_atomic_decomp(list_atomic_index,path_rfd)       
        print('   <<<< Done')
    except BaseException as e:
        
        print('Error in the Run_RF_AtomicDecomposition() method ')
        print(e)
    
    
def Run_RFC_AtomicDecomposition(path,list_atomic_index):
    
    try:
        
        path_rfcd=os.path.join(cwd,'rfcd')
        make_dir2(path_rfcd)
        
        instance=RFC_atomic_decomp(path)
        print('   <<<< Processing...')
        instance.plot_multiple_atomic_decomp(list_atomic_index,path_rfcd)
        print('   <<<< Done')
    
    except BaseException as e:
        
        print('Error in the Run_RFC_AtomicDecomposition() method ')
        print(e)
    
    
def Run_RF_Fragment_Decomposition(path,list_fragments,list_fragment_names=None):
    
    try:
        
        path_rffd=os.path.join(cwd,'rfd')
        make_dir2(path_rffd)
        
        instance=RFC_atomic_decomp(path)
        print('   <<<< Processing...')
        instance.plot_multiple_fragments_force_decomp(list_fragments,path_rffd,list_fragment_names)    
        print('   <<<< Done')
    
    except BaseException as e:
        
        print('Error in the Run_RF_Fragment_Decomposition() method ')
        print(e)
    
    
def Run_RFC_Fragment_Decomposition(path,list_fragments,list_fragment_names=None):
    
    try:
        
        path_irc=os.path.join(cwd,'IRC')
        path_rfcfd=os.path.join(cwd,'rfcd')
        make_dir2(path_rfcfd)
        
        instance=RFC_atomic_decomp(path)   
        print('   <<<< Processing...')
        instance.plot_multiple_fragments_force_constant_decomp(list_fragments,path_rfcfd,list_fragment_names)       
        print('   <<<< Done')
    except BaseException as e:
        
        print('Error in the Run_RF_Fragment_Decomposition() method ')
        print(e)
        
def Run_WBOAnalysis(cwd,ID):
    
    try:
        irc_path=os.path.join(cwd,'IRC')
        for file in os.listdir(irc_path):
            if file.endswith('log') and int(file[5:-6])==ID:
                irc_filename=os.path.join(irc_path,file)
        
        path_irc_points=os.path.join(cwd,'paths')
        key1,key2='IRCRx'+str(ID)+'p1','IRCRx'+str(ID)+'p2' 
        for key in [key1,key2]:
            if os.path.isdir(os.path.join(path_irc_points,key)):
                list_irc_points=[]
                my_folder=os.path.join(path_irc_points,key)
                for file in os.listdir(my_folder):
                    if file.endswith('log'):
                        filename=os.path.join(my_folder,file)
                        list_irc_points.append(filename)
            instance=WBO_analysis_1F(list_irc_points,irc_filename)
            output_dir=os.path.join(cwd,'wboa')
            make_dir2(output_dir)
            instance.plot(output_dir)
        
    except BaseException as e:
        
        sys.exit('Error in the RunWBOAnalysis method!!!')

    
    
def check_integer(value):
        
    try: 
        
        value=int(value)
        result=True
        
    except ValueError:        
        result=False
        
    return result

def check_IRC_validity(path):
    
    result='RELIABLE'
    
    try:
        
        instance=RFA_1file(path)
        Rc=instance.get_Rx_coords()
        RF=instance.get_Rx_force()[0]
        RE=instance.get_Energy()[0]
        beta_pos=Rc.index(0)
        Rc_part1=0
        Rc_part2=0
        for i in Rc:
            if i<0:
                Rc_part1+=1
            if i>0:
                Rc_part2+=1
        
        if RE[0]==RE[-1]:
            result='UNRELIABLE'
        if Rc_part1==Rc_part2:
            if RF[beta_pos]<-0.5 or RF[beta_pos]>0.5:
                result='UNRELIABLE'
            if RE[beta_pos]!=max(RE):
                result='UNRELIABLE'
        if Rc_part1!=Rc_part2:
            if 0 < abs(RE[0]-RE[-1])< 2:
                result='UNRELIABLE'
            if (RF[beta_pos]<-0.5 or RF[beta_pos]>0.5) and RE[beta_pos]!=max(RE):
                result='UNRELIABLE'
        for i in range(5):
            if RE[i]==max(RE):
                result='UNRELIABLE'
                break
            
    except BaseException as e:        
        result='ERROR'
        
    return result       
    
def make_dir(path):
    
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.mkdir(path)
    
def make_dir2(path):
    
    if os.path.isdir(path)==False:
        os.mkdir(path)
    
if __name__=='__main__':
    
    main()
    


# In[ ]:





# In[ ]:




