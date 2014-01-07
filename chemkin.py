#!/usr/bin/env python
# encoding: utf-8

"""
This module contains several helper functions useful when running Chemkin jobs
from the command line.
"""

import os.path
import subprocess
import csv
import math
import numpy
import logging
import sys

################################################################################

# The directory in which Chemkin is installed
CHEMKIN_DIR = os.path.abspath('/home/reaction/chemkin15131_linuxx8664')

# The preamble to each Chemkin execution shell script
CHEMKIN_SCRIPT_PREAMBLE = """#!/bin/sh -v

# Define Chemkin running environment
. {0}

""".format(os.path.join(CHEMKIN_DIR, 'bin', 'chemkinpro_setup.ksh'))

################################################################################

# Set logging parameters

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Create console handler and set level to debug; send everything to stdout
# rather than stderr
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)

logging.addLevelName(logging.CRITICAL, 'Critical: ')
logging.addLevelName(logging.ERROR, 'Error: ')
logging.addLevelName(logging.WARNING, 'Warning: ')
logging.addLevelName(logging.INFO, '')

formatter = logging.Formatter('%(levelname)s%(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

################################################################################

class ChemkinJob(object):
    """
    A single (gas-phase) Chemkin job.
    """
    
    def __init__(self, name, chemFile, tempDir):
        self.name = name
        self.chemFile = chemFile
        self.tempDir = tempDir
        
    @property
    def ascFile(self):
        return os.path.join(self.tempDir, '{0}_gas.asc'.format(self.name))
    
    @property
    def inputFile(self):
        return os.path.join(self.tempDir, '{0}.inp'.format(self.name))
    
    @property
    def outputFile(self):
        return os.path.join(self.tempDir, '{0}.out'.format(self.name))
    
    @property
    def monFile(self):
        return os.path.join('{0}.out.mon'.format(self.name))
    
    @property
    def ckcsvFile(self):
        return os.path.join(self.tempDir, 'CKSoln.ckcsv')
    
    @property
    def dataZipFile(self):
        return os.path.join(self.tempDir, 'XMLdata_{0}.zip'.format(self.name))    
    
    def preprocess(self):
        """
        Run the Chemkin preprocessor on the chemistry input file.
        """
        # Write the preprocessor execution script to a file
        scriptFile = os.path.join(self.tempDir, 'RUNJOB_CKPreProcess.sh')
        with open(scriptFile, 'w') as stream:
            stream.write(CHEMKIN_SCRIPT_PREAMBLE)
            stream.write('# Run the Chemkin preprocessor\n')
            stream.write('{0} -i {1} -o {2}_gas.out -c {2}_gas.asc\n\n'.format(
                os.path.join(CHEMKIN_DIR, 'bin', 'chem'),
                os.path.abspath(self.chemFile),
                os.path.join(self.tempDir, self.name),
            ))
            
        # Execute the preprocessor script
        process = subprocess.Popen(('/bin/sh', scriptFile), cwd=self.tempDir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            print stdout
            print stderr
            quit()
        
    def run(self, jobParams, model, pro=True):
        """
        Run a Chemkin job.
        """
        dtdFileSrc = os.path.join(CHEMKIN_DIR, 'data', 'chemkindata.dtd')
        dtdFileDest = os.path.join(self.tempDir, 'chemkindata.dtd')
        scriptFile = os.path.join(self.tempDir, 'RUNJOB_CKRunProcessor.sh')
        # Write the job parameters to a file
        with open(self.inputFile, 'w') as stream:
            stream.write(jobParams)
        # Write the job execution script to a file
        with open(scriptFile, 'w') as stream:
            stream.write(CHEMKIN_SCRIPT_PREAMBLE)
            stream.write('# Delete any intermediate files left over from a previous run\n')
            stream.write('rm -f {0}\n\n'.format(self.dataZipFile))
            stream.write('rm -f {0}\n\n'.format(dtdFileDest))
            stream.write('# Copy files to working directory\n')
            stream.write('ln -s {0} {1}\n'.format(dtdFileSrc, dtdFileDest))
            stream.write('# Run the job\n')
            if pro:
                stream.write('CHEMKIN_MODE=Pro\n')
                stream.write('export CHEMKIN_MODE\n')
                stream.write('{0} -i {1} -o {2} -x {3} Pro -m {4} -c {5}\n'.format(model, self.inputFile, self.outputFile, self.dataZipFile, self.monFile, self.ascFile))
            else:
                stream.write('{0} -i {1} -o {2} -x {3} -m {4} -c {5}\n'.format(model, self.inputFile, self.outputFile, self.dataZipFile, self.monFile, self.ascFile))
        # Execute the job script
        subprocess.call(('/bin/sh', scriptFile), cwd=self.tempDir)

    def postprocess(self, sens=False, rop=False, all=False, transpose=True, preferenceFile = False):
        """
        Run the Chemkin postprocessor.
        """
        scriptFile = os.path.join(self.tempDir, 'RUNJOB_CKPostProcessor.sh')
        # Write the postprocessor execution script to a file
        with open(scriptFile, 'w') as stream:
            stream.write(CHEMKIN_SCRIPT_PREAMBLE)
            stream.write('# Extract solution data to CKCSV\n')
            stream.write('GetSolution {1}{2}{3}{4}{0}\n\n'.format(
                self.dataZipFile, 
                '-nosen ' if not sens else '',
                '-norop ' if not rop else '',
                '-all ' if all else '',
                '-p'+' CKSolnList.txt ' if preferenceFile else '',
            ))
            if transpose:
                stream.write('# Transpose the data to CSV\n')
                stream.write('CKSolnTranspose -column 500 {0}\n\n'.format(self.ckcsvFile))
            stream.write('# Delete postprocessor log file\n')
            stream.write('rm -f {0}\n\n'.format(os.path.join(self.tempDir, 'ckpp_*.log')))
        
        # Execute the postprocessor script
        subprocess.call(('/bin/sh', scriptFile), cwd=self.tempDir)
        
 
################################################################################
        
    def writeinputRCM(self,problemType, reactants, temperature, pressure, endTime, 
                      Continuations=False, typeContinuation = None, Tlist = [], Plist = [],
                      variableVolume=False, variableVolumeProfile = None, 
                      solverTimeStepProfile = None):
        """
        Write input file for typical RCM
        """
        
        # Problem definition block 
                    
        input_stream=("""
! 
! problem type definition
!""")

        if problemType.lower() == 'constrainVandSolveE'.lower():
            input_stream+=("""
CONV   ! Constrain Volume And Solve Energy Equation
ENRG   ! Solve Gas Energy Equation""")
        elif problemType.lower() == 'constrainPandSolveE'.lower():
            input_stream+=("""
ENRG   ! Solve Gas Energy Equation""")
        elif problemType.lower() == 'constrainPandT'.lower():
            input_stream+=("""
TGIV   ! Constrain Pressure And Temperature
""")
        elif problemType.lower() == 'constrainVandT'.lower():
            input_stream+=("""
COTV   ! Constrain Volume And Temperature
TGIV   ! Fix Gas Temperature""")
             
        # Solver type definition block 
        
        input_stream+=("""
TRAN   ! Transient Solver""")
        
        # Physical property
        
        input_stream+=("""
! 
! physical property
! 
!Surface_Temperature   ! Surface Temperature Same as Gas Temperature
IFAC 0.1   ! Ignition Noise Filtering Factor
""")
            
        input_stream+=('PRES {0:g}   ! Pressure (atm)\n'.format(pressure/1.01325))
        input_stream+=('QLOS 0.0   ! Heat Loss (cal/sec)\n')
        input_stream+=('TEMP {0:g}   ! Temperature (K)'.format(temperature))
        
        if variableVolume:
            with open(variableVolumeProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    time = float(row[0].split(',')[0]) # (sec)
                    vol = float(row[0].split(',')[1]) # (cm3)
                    input_stream+=("""
VPRO {0:g} {1:g}   ! Volume (cm3)""".format(time,vol))
                       
        
        # Species property block
        
        input_stream+=("""
!
! species property
!
""")
        
        for reac , conc in reactants:            
            input_stream+=('REAC {0} {1:g} ! Reactant Fraction (mole fraction) \n'.format(reac,conc))
            
        # Solver control block
        
        input_stream+=("""! 
! solver control
! 
ADAP   ! Save Additional Adaptive Points
ASTEPS 20   ! Use Solver Integration Steps
ATLS 1.0E-6   ! Sensitivity Absolute Tolerance
ATOL 1.0E-20   ! Absolute Tolerance
MAXIT 4   ! Maximum Number of Iterations
RTLS 0.0001   ! Sensitivity Relative Tolerance
RTOL 1.0E-8   ! Relative Tolerance""")
        
        if solverTimeStepProfile:
            with open(solverTimeStepProfile, 'rb') as csvfile2:
                timeStepReader = csv.reader(csvfile2, delimiter=' ', quotechar='|')
                for row in timeStepReader:
                    time = float(row[0].split(',')[0]) # (sec)
                    vol = float(row[0].split(',')[1]) # (sec)
                    input_stream+=("""
STPTPRO {0:g} {1:g}           ! Solver Maximum Step Time (sec)""".format(time,vol))
            
        
        input_stream+=("""
TIME {0:g}                 ! End Time (sec)""".format(endTime))
        
        input_stream+=("""
! 
! output control and other misc. property
! 
EPSR 0.01   ! Threshold for Rate of Production
EPSS 0.001   ! Threshold for Species Sensitivity
EPST 0.001   ! Threshold for Temperature Sensitivity
GFAC 1.0   ! Gas Reaction Rate Multiplier
PRNT 1   ! Print Level Control
SIZE 10000000   ! Solution Data Block Size (bytes)""")
              
        if Continuations:
            if numpy.array(Tlist).size:                
                for i in range(numpy.array(Tlist).shape[0]):
                    input_stream+=("""
{0}
END
TEMP {1:g}""".format(typeContinuation,numpy.array(Tlist)[i]))
            
            if numpy.array(Plist).size:         
                for i in range(numpy.array(Plist).shape[0]):
                    input_stream+=("""
{0}
END
PRES {1:g}""".format(typeContinuation,numpy.array(Plist)[i]/1.01325))
                                            
        input_stream+=('\nEND')
        
        return input_stream      

################################################################################
    def writeinputPlugFlow(self,problemType, reactants, flowrate, 
                           startingAxialPosition, endingAxialPosition, diameter,
                           temperature  = None, pressure  = None,  
                           temperatureProfile = None, pressureProfile = None):
    
        """
        Write input file for typical plug flow
        """
        
        # Problem definition block 
                    
        input_stream=("""
! 
! problem type definition
! 
MOMEN ON   ! Turn on Momentum Equation
PLUG   ! Plug Flow Reactor
RTIME ON   ! Turn on Residence Time Calculation
""")

        if problemType.lower() == 'FixGasTemperature'.lower():
            input_stream+=("""
TGIV   ! Fix Gas Temperature""")
        elif problemType.lower() == 'solveGasEnergy'.lower():
            input_stream+=("""
ENRG   ! Solve Gas Energy Equation""")
    
        
        # Physical property
        
        input_stream+=("""
! 
! physical property
! 
!Surface_Temperature   ! Surface Temperature Same as Gas Temperature
""")
        
        input_stream+=('FLRT {0:g}   ! Mass Flow Rate (g/sec)\n'.format(flowrate))
        
        if pressureProfile:
            with open(pressureProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    pos = float(row[0].split(',')[0]) # (cm)
                    pres = float(row[0].split(',')[1])*1.0/1.01325 # (atm)
                    input_stream+=("""
PPRO {0:g} {1:g}   ! Pressure (atm)""".format(time,vol))
        else:
            input_stream+=('PRES {0:g}   ! Pressure (atm)\n'.format(pressure/1.01325))
            
        if temperatureProfile:
            with open(temperatureProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    pos = float(row[0].split(',')[0]) # (cm)
                    temp = float(row[0].split(',')[1]) # (K)
                    input_stream+=("""
TPRO {0:g} {1:g}   ! Temperature (K)""".format(pos,temp))
        else:
            input_stream+=('TEMP {0:g}   ! Temperature (K)'.format(temperature)) 
        
        # reactor dimension definition
        
        input_stream+=("""
VIS 0.0   ! Mixture Viscosity at Inlet (g/cm-sec)
! 
! reactor dimension definition
! 
""")
        input_stream+=('DIAM {0:g}   ! Diameter (cm)\n'.format(diameter))
        input_stream+=('XEND {0:g}   ! Ending Axial Position (cm)\n'.format(endingAxialPosition))
        input_stream+=('XSTR {0:g}   ! Starting Axial Position (cm)\n'.format(startingAxialPosition))
                       
        
        # Species property block
        
        input_stream+=("""
!
! species property
!
""")
        
        for reac , conc in reactants:            
            input_stream+=('REAC {0} {1:g} ! Reactant Fraction (mole fraction) \n'.format(reac,conc))
            
        # Solver control block
        
        input_stream+=("""! 
! 
! solver control
! 
ACHG 0.0   ! Maximum Absolute Change in Site Fractions
ADAP   ! Save Additional Adaptive Points
ATLS 1.0E-6   ! Sensitivity Absolute Tolerance
ATOL 1.0E-9   ! Absolute Tolerance
MAXIT 4   ! Maximum Number of Iterations
NNEG   ! Force Non-negative Solution
PSV 1.0E-8   ! Scaling Factor for Relaxing Surface Equations (cm/sec)
RCHG 1.0E-6   ! Maximum Relative Change in Site Fractions
RTLS 0.0001   ! Sensitivity Relative Tolerance
RTOL 1.0E-6   ! Relative Tolerance
TSTP 1.0   ! Initial Integration Step (cm""")          
        
        input_stream+=("""
! 
! output control and other misc. property
! 
GFAC 1.0   ! Gas Reaction Rate Multiplier""")
                                            
        input_stream+=('\nEND')
        
        return input_stream      
    



################################################################################

def getIgnitionDelay(ckcsvFile, tol=1.0, species = []):
    """
    Return the ignition delay time from the given CKCSV file. This function
    uses (dP/dt)_max to find the ignition delay time. A ValueError is raised if
    this dP/dt value is below a certain threshold.
    
    Alternatively, provide a list of species which are maximized at tau.  
    If the list contains a single species label, ie. ['OH'] the ignition delay will be given as the
    time at which the concentraiton of that species is maximized.  If a list of species is given,
    then the multiplicative max of the concentrations of those species is returned.  ie. if 
    it is desired to calculated the ignition delay using max([CH][O]), then the species list 
    ['CH','O'] should be inputted.
    """
    
    tdata = None; Pdata = None; specdata = []
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            if label.startswith('Time'):
                tdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
            elif label.startswith('Pressure'):
                Pdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                Punits = row[1].strip()[1:-1].lower()
                Pdata *= {'dyne/cm2': 0.1, 'atm': 101325., 'Pa': 1.0, 'bar': 1e5, 'torr': 133.322368, 'mtorr': 0.133322368, 'psi': 6894.75729}[Punits] * 1e-5
            for spec in species:
                specdata.append(numpy.array([float(r) for r in row[2:]], numpy.float))
    
    if tdata is None or Pdata is None:
        raise Exception('Unable to read time and/or pressure data from the given CKCSV file.')
    
    if species:
        multdata = numpy.ones(len(specdata[0]))
        for spec in specdata:
            multdata *= spec
        
        index = multdata.argmax()
        return tdata[index]
    
    else:
        dPdt = (Pdata[1:] - Pdata[:-1]) / (tdata[1:] - tdata[:-1])
        dPdt = dPdt[numpy.isfinite(dPdt)]
         	
        index = dPdt.argmax()
        if dPdt[index] < tol:
            raise ValueError('No ignition occurred in the given simulation.')
        
        return 0.5 * (tdata[index] + tdata[index+1])

################################################################################

def runIgnitionReactionSensitivity(runChemkinJob, inputFile, dictionaryFile):
    """
    Supply a runChemkinJob python function which returns the ignition delay with a 
    chemkin file input.  This will run finite difference reaction sensitivities and 
    save them to a csv file.
    """
    
    from rmgpy.chemkin import loadChemkinFile, saveChemkinFile
    
    
    speciesList, reactionList = loadChemkinFile(inputFile, dictionaryPath = dictionaryFile, readComments = False)
    
    num_reactions = len(reactionList)
    
    factor_high = 1.05
    factor_low = 1. / factor_high
    
    worksheet = csv.writer(file('ignition_rxn_sensitivity.csv', 'w'))
    worksheet.writerow(['Index', 'Reaction', 'd[ln k]','tau_high','tau_low','d[ln tau]/d[ln k]'])
    
    logging.info('Running reaction sensitivity analysis using finite differences...')
    for index, reaction in enumerate(reactionList):
        rxn_index = index + 1
        rxn_string = reaction.toChemkin(kinetics = False)
        logging.info('At reaction {0} of {1}. {2}'.format(rxn_index, num_reactions, rxn_string))
                
        reaction.kinetics.changeRate(factor_high)
        saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose = False)
        tau_high = runChemkinJob('chem_temp.inp')
        reaction.kinetics.changeRate(1./factor_high)    # reset the kinetics
        
        reaction.kinetics.changeRate(factor_low)
        saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose = False)
        tau_low = runChemkinJob('chem_temp.inp')
        reaction.kinetics.changeRate(1./factor_low)     # reset the kinetics
        
        if tau_high != 0 and tau_low != 0:
            sens = numpy.log(tau_high / tau_low) / numpy.log(factor_high / factor_low)
        else:
            sens = 0
            
        worksheet.writerow([rxn_index, rxn_string, numpy.log(factor_high / factor_low), tau_high, tau_low, sens])
        
################################################################################

def runIgnitionThermoSensitivity(runChemkinJob, inputFile, dictionaryFile):
    """
    Supply a runChemkinJob python function which returns the ignition delay with a 
    chemkin file input.  This will run finite difference sensitivities to enthalpies and 
    save them to a csv file.
    """
    
    from rmgpy.chemkin import loadChemkinFile, saveChemkinFile, getSpeciesIdentifier
    from rmgpy.quantity import Quantity
    
    
    speciesList, reactionList = loadChemkinFile(inputFile, dictionaryPath = dictionaryFile, readComments = False)
    
    num_species = len(speciesList)
    
    deltaH = Quantity(0.5, 'kcal/mol').value_si
    
    worksheet = csv.writer(file('ignition_thermo_sensitivity.csv', 'w'))
    worksheet.writerow(['Species', 'd[del H] (kcal/mol)', 'tau_high', 'tau_low', 'd[ln tau]/d[del H]'])
    
    logging.info('Running thermo sensitivity analysis using finite differences...')
    for index, species in enumerate(speciesList):
        species_index = index + 1
        species_string = getSpeciesIdentifier(species)
        logging.info('At species {0} of {1}. {2}'.format(species_index, num_species, species_string))
        
        species.thermo.changeBaseEnthalpy(deltaH)
        saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose = False)
        tau_high = runChemkinJob('chem_temp.inp')
        species.thermo.changeBaseEnthalpy(-deltaH)    # reset the thermo
        
        species.thermo.changeBaseEnthalpy(-deltaH)
        saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose = False)
        tau_low = runChemkinJob('chem_temp.inp')
        species.thermo.changeBaseEnthalpy(deltaH)     # reset the kinetics
        
        if tau_high != 0 and tau_low != 0:
            sens = numpy.log(tau_high / tau_low) / (2 * deltaH)
        else:
            sens = 0
            
        worksheet.writerow([species_string, '1', tau_high, tau_low, sens])
      
