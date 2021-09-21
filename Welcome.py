import time

def message():
    
    time.sleep(1)
    name='AMADAR : Automated workflow for Mechanistic Analysis of the Diels-Alder Reaction'
    size=len(name)
    print('+-'*(13+size//2))
    title1='\nCOMPUTATIONAL MECHANISTIC CHEMISTRY AND DRUG DISCOVERY'
    print(title1+'\n')
    print('Rhodes University')
    print('Department of Chemistry')    
    print(name+'\n')
    print('+-'*(13+size//2))
    time.sleep(1)
    print('\nThis code is mainly designed to:')
    print('     1. Automate the generation of Diels-Alder transition states geometries')
    print('     2. Determine IRC paths from the predicted TS(s) following the details in the da.ini file')
    print('     3. Perfom several analyses based on the IRC paths, following the details in the analysis.ini file')
    print('     4. Optimize reactants and cycloadducts')
    print('     5. Carry out the retro-DA transformation of the cycloadducts')
    print('\nPlease report any issue to :')
    print('     k.lobb@ru.ac.za')
    print('     isamurabft@gmail.com')
    print('#'*107)
    print(time.ctime())
    print('\n')    
    time.sleep(2)
    msg='PROGRESS REPORT'
    s=len(msg)
    print(msg)
    print('+-'*(1+s//2))
    print('\n')