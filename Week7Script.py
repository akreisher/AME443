
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

from scipy import signal

import matplotlib.pyplot as plt
import numpy as np
import os, csv


print('\n Week 7 Script')
print('-----------')

week3Files = (f for f in os.listdir('Week 3') if f[-3:] == 'csv')
freeRuns = (FreeResponse(f,3) for f in week3Files if f.split('_')[0]=='Free')
free484 = next(run for run in freeRuns if run.mass == .484)
free969 = next(run for run in freeRuns if run.mass == .969)
free484.SystemID()
free969.SystemID()

free484.parameters(free969)

print(free484.zeta)

stepCLFile = 'c:/Users/Alex/Documents/AME443/Week 7/Step_m484_k3k3_r100_cl_1.csv'
stepCL = StepResponse(stepCLFile)
stepCL.SystemID()

mTot = free484.m+stepCL.mass

KHWCL = free484.k*stepCL.xInf/stepCL.stepInput

zeta = free484.zeta
omegaN = free484.omegaN

omegaMin = 1.5*omegaN
omegaMax = 2*omegaN
zetaMax = 0.2;
zetaMin = 0.5;

KpMin = (omegaMin**2-omegaN**2)/(KHWCL/mTot);
#KpMax = (omegaMax**2-omegaN**2)/(KHWCL/mTot);
KpMax = 0.35;

KdRange = np.linspace(0,0.01,10000)
KdMin = (4*(omegaN**2)*((zeta**2)-zetaMax**2)+4*zeta*omegaN*(KHWCL/mTot)*KdRange + (KHWCL*KdRange/mTot)**2)/(4*KHWCL*(zetaMax**2)/mTot);
KdMax = (4*(omegaN**2)*((zeta**2)-zetaMin**2)+4*zeta*omegaN*(KHWCL/mTot)*KdRange + (KHWCL*KdRange/mTot)**2)/(4*KHWCL*(zetaMin**2)/mTot);


plt.figure(1)
plt.plot([0,1],[KpMin,KpMin],label='Kp Min')
plt.plot([0,1],[KpMax,KpMax],label='Kp Max')
plt.plot(KdRange,KdMin,label='Kd Min')
plt.plot(KdRange,KdMax,label='Kd Max')
plt.xlabel('Kd')
plt.ylabel('Kp')
plt.xlim([0,0.0075])
plt.ylim([0,0.5])
plt.legend()
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

goodList = [[0.0018,0.1691],[.005,.1691],[.0024,.35],[.0063,.35]]

for ind,point in enumerate(goodList):
   Kd=point[0]
   Kp=point[1]
   path = 'Step_m484_k3k3_r500_PD_KP'+str(Kp)+'_KD'+str(Kd)+'_1.csv'
   raw = np.genfromtxt('Week 7/good/'+path,delimiter=',',skip_header=1)
   t = raw[:,0]
   t = t[t<5]
   pos = raw[:len(t),9]


   omegaNCL = np.sqrt(omegaN**2 + Kp*KHWCL/mTot)
   zetaCL = (2*zeta*omegaN + KHWCL*Kd/mTot)/(2*omegaNCL)

   num = Kp*KHWCL/mTot
   den = [1,2*zetaCL*omegaNCL,omegaNCL**2]

   tf = signal.TransferFunction(np.array(num),np.array(den))
   tsim,ysim = signal.step(tf,T=t)

   
   plt.figure(ind+2)
   plt.title('Kp = '+str(Kp)+' Kd = '+str(Kd))
   plt.plot(t,pos,label='data')
   plt.plot(tsim,ysim*500,'--',label='model')
   plt.legend()
   plt.xlim([0,0.6])
   ax = plt.gca()
   ax.spines['top'].set_visible(False)
   ax.spines['right'].set_visible(False)

figs = [plt.figure(n) for n in plt.get_fignums()]

for ind,fig in enumerate(figs):
    fig.savefig('Week 7/Checkin/Fig'+str(ind+1)+'.png')
plt.show()









