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
path = "Week7/good/"
files = os.listdir(path)

for file in files:
    print(file)
runs = (StepResponse(path+file,5) for file in files if file.split('_')[0]=="Step")
for run in runs:
    
    run.plot()
    run.SystemID(method = 'overshoot')
    print(file)
    print("ts = ", run.tSettling)
    print("overhsoot = ", run.overshoot)
    print()    
