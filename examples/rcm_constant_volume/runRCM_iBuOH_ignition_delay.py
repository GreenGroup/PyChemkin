#!/usr/bin/env python
# encoding: utf-8

"""
This script uses Chemkin to run simulations of a rapid compression machine 
using an RMG-generated model. 
"""

import os.path
import numpy
import pylab
import matplotlib.ticker

from chemkin import ChemkinJob, getIgnitionDelay

currentDir = os.path.dirname(__file__)

job = ChemkinJob(
    name = 'iBuOH',
    chemFile = os.path.join(currentDir, 'C4H9OH_NSRL_model.inp'),
    tempDir = os.path.join(currentDir, 'temp'),
)

# Solver time step profile
#solTimeStepFile=os.path.realpath(os.path.join(currentDir, 'RCM','time_step.csv'))

job.preprocess()

Tlist = 1000./numpy.arange(1, 1.4, 0.05)
P = 15.     # Pressure (bar)

iBuOH = numpy.array([0.0338,0.0338,0.0338])
o2 = numpy.array([0.2030,0.4060,0.1015])
n2 = numpy.array([0.7632,0.5602,0.8647])
time = 0.5

taulist = numpy.zeros((Tlist.shape[0],iBuOH.shape[0]))


for j in range(iBuOH.shape[0]):
    
    for i in range(Tlist.shape[0]):   
    
        input = job.writeInputHomogeneousBatch(
                                  
                                  problemType = 'constrainVandSolveE', # Problem type(constrainVandSolveE,
                                                                       # constrainPandSolveE, constrainPandT, constrainVandT)                                      
                                  
                                  reactants=[('N2',n2[j]),
                                             ('iC4H9OH',iBuOH[j]),
                                             ('O2',o2[j]),],     # Reactant (mole fraction)
                                  
                                  temperature = Tlist[i], # Temperature(K)
                                  pressure = P ,   # Pressure (bar)
                                  endTime = time ,   # End Time (sec)
                                  
                                  #Continuations = True,             # Continuations
                                  #typeContinuation = 'NEWRUN',      # Type of continuation NEWRUN or CNTN
                                  #Tlist=[1000.,1200.],              # Temperature (K) list of continuations
                                  #Plist=[1.,2.],                   # Pressure (atm) list of continuations
                                  
                                  #variableVolume = True,            # Variable volume true / false
                                  #variableVolumeProfile =vproFile,  # File name path for vpro (sec,cm3)
                                  
                                  #solverTimeStepProfile = solTimeStepFile # Solver time profile (sec)
                                  )
        job.run(input, model='CKReactorGenericClosed', pro=True)
        job.postprocess(sens=False, rop=False, all=True, transpose=False)
    
        try:
            taulist[i,j] = getIgnitionDelay(job.ckcsvFile)
        except ValueError:
            pass
      
    print '========= =========='
    print 'T (K)     tau (ms)  '
    print '========= =========='
    for i in range(Tlist.shape[0]):
        if taulist[i,j] != 0:
            print '{0:9g} {1:9.3g}'.format(Tlist[i], taulist[i,j]*1000.)
    print '========= =========='

data_dump =numpy.concatenate([numpy.array([Tlist]).T,taulist],axis=1)
numpy.save('data_notair_15bar.npy',data_dump)

fig = pylab.figure(figsize=(5,4))
pylab.semilogy(1000./Tlist, taulist[:,0] * 1000., '-b')
pylab.semilogy(1000./Tlist, taulist[:,1] * 1000., '-g')
pylab.semilogy(1000./Tlist, taulist[:,2] * 1000., '-r')

#pylab.xlim(0.9999, 1.8001)
pylab.xlabel('1000 / (Temperature (K))')
pylab.ylabel('Ignition delay (ms)')

ax = pylab.gca()
ax.xaxis.set_major_locator( matplotlib.ticker.MultipleLocator(0.1) )
ax.xaxis.set_minor_locator( matplotlib.ticker.MultipleLocator(0.05) )

fig.subplots_adjust(left=0.12, bottom=0.11, right=0.95, top=0.94, wspace=0.33, hspace=0.35)
pylab.savefig('ignition_delays_RCM.pdf')
pylab.savefig('ignition_delays_RCM.png')




