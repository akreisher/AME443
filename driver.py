from importlib import reload
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
#freeList = [FreeResponse(file,3) for file in os.listdir('Week 3') if file.split('_')[0]=='Free' and file[-3:]!='pdf']
stepList = [StepResponse(file,3,20) for file in os.listdir('Week 3') if file.split('_')[0]=='Step' and file[-3:]!='pdf']
sweepList = [SineSweep(file,4) for file in os.listdir('Week 4') if file.split('_')[0]=='SineSweep' and file[-3:]!='pdf']
constList = [ConstSine(file,4) for file in os.listdir('Week 4') if file.split('_')[0]=='Const']

#free484 = [obj for obj in freeList if int(obj.mass*1000)==484]

#free969 = [obj for obj in freeList if obj not in free484]


step484 = [obj for obj in stepList if int(obj.mass*1000)==484]
step969 = [obj for obj in stepList if obj not in step484]
sweep484 = [obj for obj in sweepList if int(obj.mass*1000)==484]
sweep969 = [obj for obj in sweepList if obj not in sweep484]



n=len(sweep484)

zetas = np.zeros(n,)
omegaNs = np.zeros(n,)
ms = np.zeros(n,)
cs = np.zeros(n,)
ks = np.zeros(n,)
i=0
for run484,run969 in zip(sweep484,sweep969):
	out = run484
	other = run969
	method = 'hp'
	run484.SystemID(method)
	run969.SystemID(method)
	out.parameters(other)
	zetas[i]=out.zeta
	omegaNs[i]=out.omegaN
	ms[i] = out.m
	cs[i]=out.c
	ks[i] = out.k
	i+=1
zeta = np.mean(zetas)
omegaN = np.mean(omegaNs)
m = np.mean(ms)
c = np.mean(cs)
k = np.mean(ks)
print("zetas = ",zetas)
print("zeta = ",zeta)
print("\n")
print("omegaNs = ", omegaNs)
print("omegaN = ", omegaN)
print("\n")
print("ms = ", ms)
print("m = ", m)
print("\n")
print("cs = ", cs)
print("c = ", c)
print("\n")
print("ks = ", ks)
print("k = ", k)



"""
freeZetas = np.empty([len(freeList),1])
freeOmegas = np.empty(np.shape(freeZetas))
for i,obj in enumerate(freeList):
	obj.SystemID('log')
	freeZetas[i] = obj.zeta
	freeOmegas[i] = obj.omegaN
print("free zetas = ",freeZetas)
print("free omegas = ",freeOmegas)


stepZetas = np.empty([len(stepList),1])
stepOmegas = np.empty(np.shape(stepZetas))
for i,obj in enumerate(stepList):
	obj.SystemID('rise')
	stepZetas[i] = obj.zeta
	stepOmegas[i] = obj.omegaN
	plotargs = {"fig" : 1, "data" : True, "models" : 'all', "show":True};
	#obj.plot(**plotargs)
print("step zetas = ",stepZetas)
print("step omegas = ",stepOmegas)

Week4Files = [file for file in os.listdir("Week 4") if file[-3:]!='pdf']

sweepList = []
constList = []
for file in Week4Files:
	params = file.split("_")
	if params[0] == 'SineSweep':
		sweepList.append(SineSweep(file,4))
	elif params[0] == 'Const':
		constList.append(ConstSine(file,4))


sweepZetas = np.empty([len(sweepList),1])
sweepOmegas = np.empty(np.shape(sweepZetas))

for i,obj in enumerate(sweepList):
	obj.SystemID('hp')
	sweepZetas[i] = obj.zeta
	sweepOmegas[i] = obj.omegaN
	plotargs = {"fig" : 1, "toPlot" : 'gain','models':'all', 'show': True};
	#obj.plot(**plotargs)
print("sweep zetas = ",sweepZetas)
print("sweep omegas = ",sweepOmegas)

constZetas = np.empty([len(constList),1])
constOmegas = np.empty(np.shape(constZetas))

for i,obj in enumerate(constList):
	obj.SystemID('hp')
	constZetas[i] = obj.zeta
	constOmegas[i] = obj.omegaN
	plotargs = {"fig" : 1, "toPlot" : 'bode','models':'all', 'show': True};
	obj.plot(**plotargs)
print("const zetas = ",constZetas)
print("const omegas = ",constOmegas)
"""
