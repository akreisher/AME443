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


print('\nTest Script')
print('-----------')

stepOL484File = 'c:/Users/Alex/Documents/AME443/Week 3/Step_m484_k3k3_r1.5_1.csv'
stepOL969File = 'c:/Users/Alex/Documents/AME443/Week 3/Step_m969_k3k3_r1.5_1.csv'
stepCLFile = 'c:/Users/Alex/Documents/AME443/Week 7/Step_m484_k3k3_r1000_cl_1.csv'

test = 'Step_m484_k3k3_r1000_cl_1.csv'

test.split

stepOL484 = StepResponse(stepOL484File,20)
stepOL969 = StepResponse(stepOL969File,20)
stepOL484.SystemID('fit')
stepOL969.SystemID('fit')
stepOL484.parameters(stepOL969)

m = stepOL484.m
print(m)

stepCL = StepResponse(stepCLFile)
stepCL.SystemID('fit')
print(stepCL.xInf)

#exec(open("test.py").read(), globals())


