# AMADAR

## About

ADAMAR can identify the transition state of a Diels-Alder reaction, given the product, with a 95% success rate. All that is required is an input SMILES string of the product, and bothr reactant and product 3d structures will be generated, and the transitions tate will be generated and refined at a chosen high level of theory. Further ADAMAR can analyse the IRC for the specific DA process, for example in terms of RFA.

## Getting Started

ADAMAR includes two main scripts, **--init.pynb** and **myIRCAnalyzer.ipynb**, dependent on **RDKit**, **os**,**sys**,**time** and **shutil**. There are several modules which are called by these two scripts. In addition **Gaussian** must be installed and accessible for calculation on the local machine. Support for **GAMESS** is not yet available.
There are two configuration files: **da.ini** and **analysis.ini**, in which all Diels-Alder transition state searches and analyses may be customized.


## Configuration Files

### da.ini

This file must be edited prior to generation of transition states. There are two sections to this configuration script **[job_details]** and **[flags]**. Please also see the sample da.ini provided.

**[job_details]**
NPROCSHARED = number of processors to be used for the jobs

D_SPLITTING = Final distance between the two fragments during the constrained optimization

CALC_LEVELS_RC = Levels of theory to be used in the optimization of reactants and cycloadducts. 
		Two levels must be given, separated by a semi-colon. 
		Users are advised to choose a semi-empirical and QM levels as the 1st and 2nd levels respectively.
		The defaults are PM6 and b3lyp/6-31G(d).

CALC_LEVELS_TS =Levels of theory to be used in the optimization of reactants and cycloadducts. 
		Three levels must be given, separated by a semi-colon. 
		The first level is used in the constrained optimization of the cycloadduct twoards the pseudo-guess TS
		The second level is for the refinement of the pseudo-guess into a guess TS. It is better using the same level of the 1st.
		The third level is needed for the final refinement of the guess-TS to the targeted QM accuracy.
		The defaults are PM6, PM6 and b3lyp/6-31G(d).

CALC_LEVEL_IRC = The level of theory to used in IRC calculations. This should be the same as the one used to refine the guess into the likely TS.

NBR_IRC_POINT_PP = Maximum number of points to use in the construction of the IRC path. 60 is good starting choice for mid-sized molecules, but may be increased or decreased depending on the system.

IRC_STEP_SIZE = Step-size to consider in the determination of the IRC path. 
		The default is 8, corresponding to 0.08 (amu) ^0.5 Bohr

NBR_PATHS=	This keyword can have two states : 1 in case we want to built a unique IRC file or 2 if we want to have 2 separate input files.
		The option 2 has not be activated yet, but will be very soon since the rest of analyses are still based on a unique IRC file containing both directions.

TS_ID_numbers = ID numbers of the reactions for which  for which the IRC path has to be generated. 
		These must be integers values separated by commas.
		This flag can also be set to 'ALL',which means all the TS predicted will be involved in IRC calculations.
        

**[flags]**
SCRATCH = 	This flag can only have 0 and 1 values. 
		1 tells the initialization code to regenerate the very initial input files saved in the R,C and TS folders
		0 tells the initialization code to skip the regeneration. This option will only work if there are already input files  in the R, C and TS or IRC folders.

RC_FLAG = 	This flag can only have 0 and 1 values.
		Value 1 tells the code to optimize reactants and cycloadducts, whereas 0 says the opposite.

TS_FLAG = 	This flag can only have 0 and 1 values.
		Value 1 tells the code to optimize the TSs, whereas 0 says the opposite.

IRC_FLAG = 	This flag can only have 0 and 1 values.
		Value 1 tells the code to run IRC calculations, whereas 0 says the opposite.
		This will work if and only if there are TS already predicted and saved in the ...TS/GTS/TS folder.

IRC_GEOMS_CONSTR = 	This flag can only have 0 and 1 values.
			Value 1 tells the code to generate geometries of the system from IRC paths.
			This will work only if there are IRC paths already constructed and saved in the IRC folder.
			Outputs are saved in the "paths" folder.

RFA_FLAG = 		This flag can only have 0 and 1 values.
			Value 1 is for authorizing the reaction force analysis to be done.
			Supplementary flags (keywords) should be given in the analysis.ini file to indicate which system to analyze.

RFD_FLAG = 		This flag can only have 0 and 1 values.
			Value 1 is for authorizing the atomic resolution of energy derivatives (along the IRC path) to be done.
			Supplementary flags (keywords) should be given in the analysis.ini file to indicate which system to analyze or which atoms (fragments) to invole in the decomposition.

WBOA_FLAG = 		This flag can only have 0 and 1 values.
			Value 1 enables the Wiberg Bond Order analysis to be carried out, while 0 is for the opposite.
			Supplementary flags (keywords) should be given in the analysis.ini file to indicate which system to analyze.
			It will work only if IRC natural population calculations have already been perfomed on the IRC geometries. These are saved in the "paths" folder.


OVERWRITE_FLAG = 	This flag tell whether the managder3 module to overwrite IRC folder as a new call is made, or not.
			The value 1 allows for overwriting, while 0 blocks it.
      
### analysis.ini

This file must be edited prior to analyzing IRC's that have already been generated. There are three sections to this configuration file, **[RFA]**, **[RFD]** and **[IRC_PATHS]**. Please also look at the example analysis.ini provided.

**[RFA]**
Unq_RFA = 	This keyword indicates the list of ID numbers of IRC files to include in the reaction force analysis.
		Values must be integers, separated by commas. 
		In case all the files are to be considered, the value "-1" must be used. 
		The value 0 means that no file is to be analyzed.

Multiple_RF = 	This keyword indicates the list of ID numbers of IRC files to include in the multiple RFA in order to superimpose their reaction force curves. 
		Values must be integers, separated by commas. 
		The value "-1" is not allowed to avoid that the case where a user should try to overlay 1000000000 systems :)
		The value 0 disables the analysis.

Multiple_RE = 	This keyword indicates the list of ID numbers of IRC files to include in the multiple RFA in order to superimpose their potential energy curves. 
		Values must be integers, separated by commas. 
		The value "-1" is not allowed to avoid that the case where a user should try to overlay 1000000000 systems :).
		The value 0 disables the analysis.


Multiple_RFC =  This keyword indicates the list of ID numbers of IRC files to include in the multiple RFA in order to superimpose their reaction force constant curves. 
		Values must be integers, separated by commas. 
		The value "-1" is not allowed to avoid that the case where a user should try to overlay 1000000000 systems :). 
		The value 0 disables the analysis.


**[RFD]**
JOB_ID = 	This keyword indicates the ID number of IRC file to consider for the the atomic decomposition of energy derivatives . 
		The value 0 disables the analysis.

ATOMS = 0	This keyword is for the list of atomic indexes to involve in the decomposition.
		Values must be integers, separated by commas.


FRAG =		This keyword is for the list of fragments to involve in the decomposition.
		Fragments are separated by semi-colons, while indexes of atoms from the same fragment are separated by commas.
		Make sure all the values are integers.
		The value 0 disables this analysis.

FRAG_NAMES = 	The names of the fragments must be given in the same order as indicated in the FRAG keyword.
		This section can be left blank or put any value when FRAG=0 because no fragment decomposition will be made.	


**[IRC_PATHS]**
IRC_ID= 	ID number of the IRC file for which the geometries have to be extracted and run (single point calculation,with pop analysis)

LEVEL = 	Level of theory to use for the previous job.

WBOA_ID= 	ID number of the system (IRC file) for which a Wiberg Bond order analysis should be done based on the output of the previous step.
		This will work if single point calculations for the geoemtries extracted from the IRC path of the system of interest was be run and saved in the "path" folder.
		
