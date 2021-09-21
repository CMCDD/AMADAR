import os
from configparser import ConfigParser

cwd=os.getcwd()

def RunGaussian(filename):
   
    
    config=ConfigParser()
    config.read(os.path.join(cwd,"da.ini"))
   
    os.environ['GAUSS_SCR']=config['gaussian_environment']['GAUSS_SCR']
    os.environ['GAUSS_SCRDIR']=config['gaussian_environment']['GAUSS_SCRDIR']
    os.environ['G09BASIS']=config['gaussian_environment']['G09BASIS']
    os.environ['GAUSS_ARCHDIR']=config['gaussian_environment']['GAUSS_ARCHDIR']
    os.environ['GAUSS_BSDDIR']=config['gaussian_environment']['GAUSS_BSDDIR']      
    os.environ['GAUSS_EXEDIR']=config['gaussian_environment']['GAUSS_EXEDIR']      
    os.environ['GAUSS_LEXEDIR']=config['gaussian_environment']['GAUSS_LEXEDIR']      
    os.environ['LD_LIBRARY_PATH']=config['gaussian_environment']['LD_LIBRARY_PATH']      
    os.environ["PATH"] += config['gaussian_environment']['PATH']      
   
    logname=filename[:-3]+'log'
   
    os.system('g09 < '+filename+' > '+logname)
   
    

       
