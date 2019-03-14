from importlib import reload
#Import and reload
import AME443Utility.FreeResponse
AME443Utility.FreeResponse = reload(AME443Utility.FreeResponse)
from AME443Utility.FreeResponse import FreeResponse

import AME443Utility.StepResponse
AME443Utility.StepResponse = reload(AME443Utility.StepResponse)
from AME443Utility.StepResponse import StepResponse

import AME443Utility.SineSweep
AME443Utility.SineSweep = reload(AME443Utility.SineSweep)
from AME443Utility.SineSweep import SineSweep

import AME443Utility.ConstSine
AME443Utility.ConstSine = reload(AME443Utility.ConstSine)
from AME443Utility.ConstSine import ConstSine

import matplotlib.pyplot as plt
import numpy as np
import os, csv

#Open all data files
freeList = [FreeResponse(file,3) for file in os.listdir('Week 3') if file.split('_')[0]=='Free' and file[-3:]!='pdf']
stepList = [StepResponse(file,3,20) for file in os.listdir('Week 3') if file.split('_')[0]=='Step' and file[-3:]!='pdf']
sweepList = [SineSweep(file,4) for file in os.listdir('Week 4') if file.split('_')[0]=='SineSweep' and file[-3:]!='pdf']
constList = [ConstSine(file,4) for file in os.listdir('Week 4') if file.split('_')[0]=='Const']

with open('results/FreeIDResults.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(['run','zeta','omegaN','m','k','c'])
    for run in freeList:
    	run.SystemID()
    	for other in freeList:
    		try:
                    other.SystemID()
                    run.parameters(other)
    		except ValueError:
    			continue
    		break

    	writer.writerow([run.file,run.zeta,run.omegaN,run.m,run.k,run.c])
with open('results/StepIDResults.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(['run','zeta','omegaN','m','k','c','KHW_OL'])

    for run in stepList:
    	run.SystemID()
    	for other in stepList:
    		try:
    			other.SystemID()
    			run.parameters(other)
    		except ValueError:
    			continue
    		break

    	writer.writerow([run.file,run.zeta,run.omegaN,run.m,run.k,run.c,run.K_HW])
