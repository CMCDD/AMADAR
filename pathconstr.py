import os
import sys
from IRC import IRC_1F
from configparser import ConfigParser
from Gaussian import RunGaussian
from multiprocessing import Pool
import time




def Constr_irc_paths(cwd):
    
    try:
        irc_path=os.path.join(cwd,'IRC')
        if os.path.isdir(irc_path):            
            config=ConfigParser()
            config.read(os.path.join(cwd,"analysis.ini"))
            IRC_ID=[int(x) for x in  (config['IRC_PATHS']['IRC_ID']).split(',')]
            level_of_theory=config['IRC_PATHS']['LEVEL']
            path_geoms=os.path.join(cwd,'paths')
            make_dir2(path_geoms)
            if len(IRC_ID)==1 and IRC_ID[0]==-1:
                print('\nGenerating IRC geometries for all the systems in ',irc_path)
                for file in os.listdir(irc_path):
                    if file.endswith('.log'):                        
                        filename=os.path.join(irc_path,file)
                        instance=IRC_1F(filename)
                        new_dir=os.path.join(path_geoms,file[:-4])
                        make_dir2(new_dir)
                        instance.generate_irc_geom_spcalc_com_files(new_dir,level_of_theory)
                        print(filename,' DONE')
                        time.sleep(0.50)   
                        
                        myinput_files=[] 
                        for inputfile in os.listdir(new_dir):
                            filepath=os.path.join(new_dir,inputfile)
                            myinput_files.append(filepath)
                            
                        pool=Pool(processes=4)
                        result=pool.map(RunGaussian,myinput_files)
                        pool.close()
                        pool.join()

            else:
                print('\nGenerating IRC geometries of selected systems as defined in the analysis.ini file')
                for file in os.listdir(irc_path):
                    if file.endswith('.log') and int(file[5:-6]) in IRC_ID:
                        filename=os.path.join(irc_path,file)
                        instance=IRC_1F(filename)
                        new_dir=os.path.join(path_geoms,file[:-4])
                        make_dir2(new_dir)
                        instance.generate_irc_geom_spcalc_com_files(new_dir,level_of_theory)
                        print(filename,' DONE')
                        time.sleep(0.25)  
                        
                        print('\n   <<<< Files generated')
                        myinput_files=[] 
                        for inputfile in os.listdir(new_dir):
                            filepath=os.path.join(new_dir,inputfile)
                            print(filepath)
                            time.sleep(0.10)
                            myinput_files.append(filepath)
                        
                        
                        pool=Pool(processes=4)
                        result=pool.map(RunGaussian,myinput_files)
                        pool.close()
                        pool.join()
        else:
            sys.exit('Aborted!!! Could not find IRC folder in the root directory')
            
    
    except BaseException as e:
        
        print(e)       
        sys.exit('Unable to generate IRC geometries')
        
    
def make_dir(path):
    
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.mkdir(path)
    
def make_dir2(path):
    
    if os.path.isdir(path)==False:
        os.mkdir(path)    
    
    
