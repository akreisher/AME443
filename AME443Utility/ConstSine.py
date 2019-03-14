import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit, fsolve
import os, sys, math

class ConstSine:
	def __init__(self,folder, weekNum):
		self.mass = float(folder.split("_")[-1])/1000
		dirs = os.listdir('Week '+str(weekNum));
		self.folder='Week '+str(weekNum)+'/'+folder
		self.datasets = os.listdir(self.folder)
		self.omegas = np.empty([len(self.datasets),1])
		self.dB = np.empty(np.shape(self.omegas))
		self.phase = np.empty(np.shape(self.omegas))
		self.dt = np.empty(np.shape(self.omegas))
		for i,file in enumerate(self.datasets):
			f = float(file[1:-4])
			self.omegas[i] = f*2*math.pi
			data = np.genfromtxt(self.folder+'/'+file, delimiter=",",skip_header=1)
			t = data[:,0]
			pos = data[:,9]
			ce = data[:,4]
			pksPos,_ = find_peaks(pos)
			pksCE,_ = find_peaks(ce)

			dt = t[pksPos[-1]]-t[pksCE[-1]]

			"""
			if dt<0:
				plt.plot(t,pos)
				plt.scatter(t[pksPos],pos[pksPos])
				plt.plot(t,ce)
				plt.scatter(t[pksCE],ce[pksCE])
				plt.show()
			"""

			self.dt[i] = dt
			self.phase[i] = -dt*self.omegas[i]
			self.dB[i] = 20*np.log10(pos[pksPos[-1]]/0.5)

			
		ind = np.argsort(self.omegas,axis=0)
		self.omegas = np.squeeze(self.omegas[ind])
		self.dB = np.squeeze(self.dB[ind])
		self.phase = np.squeeze(self.phase[ind])*180/math.pi
		self.zeta = None
		self.m = None
		
	def SystemID(self,method='max', retLabel = False):
		
		dB=self.dB
		omegas = self.omegas
		dB = dB-dB[0]
		maxdB = dB[0]
		indMin = None
		indMax = 0
		indPlus = None
		i=0

		label = None

		for power in dB:
			if power > maxdB:
				maxdB=dB[i]
				indMax = i
			i+=1
		self.omegaN = omegas[indMax]
		r = omegas/omegas[indMax]
		if method == 'hp':
			i=0
			for power in dB:
				if indMin is None and dB[i]>dB[indMax]-3:
					indMin = i
					rMin = r[i-1]+(r[i]-r[i-1])*(dB[indMax]-3-dB[i-1])/(dB[i]-dB[i-1])
				if indMin is not None and dB[i]<dB[indMax]-3:
					indPlus = i
					rMax = r[i-1]+(r[i]-r[i-1])*(dB[indMax]-3-dB[i-1])/(dB[i]-dB[i-1])
					break
				i+=1
			self.zeta = (rMax-rMin)/2
			label = 'Half-power'
		elif method =='max':
			self.zeta = 1/(2*(10**(maxdB/20)))
			label = 'Max magnification'
		else:
			print('System ID method must be hp or max')
		if retLabel:
			return label


	def plot(self,fig=None,show=True, toPlot='bode',models=None, col=None):
		
		omegas = self.omegas
		dB=self.dB
		
		if toPlot =='bode':
			#Plot gain
			f, (ax1,ax2) = plt.subplots(2,1,num=fig)
			ax1.scatter(self.omegas,self.dB,label='Experimental data',s=12)
			ax1.set_xscale('log')
			ax1.spines['top'].set_visible(False)
			ax1.spines['right'].set_visible(False)
			ax1.set_ylabel("Gain [dB]")
			if models is not None:
				if models == 'both' or models=='all':
					models = 'max hp'
				models = models.split()
				for model in models:
					label = self.SystemID(method=model,retLabel=True)
					gain = 10**(self.dB[0]/20)*self.omegaN**2/np.sqrt((self.omegaN**2-omegas**2)**2+(2*self.zeta*self.omegaN*omegas)**2)
					ax1.semilogx(omegas,20*np.log10(gain),label=label+' model')
					phase = (180/math.pi)*np.arctan2(-2*self.zeta*self.omegaN*omegas,self.omegaN**2-omegas**2)
					ax2.semilogx(omegas,phase,label=label+' model')
			#Plot phase
			ax2.scatter(self.omegas,self.phase,s=12)
			ax2.set_xscale('log')
			ax2.set_ylabel('Phase [deg]')
			ax2.set_ylim([-180,20])
			ax2.set_yticks(np.arange(-180,30,30))

			plt.xlabel("Frequency [rad/s]")
			ax1.legend()



		elif toPlot=='gain':
			plt.scatter(self.omegas,self.dB,label = 'Constant Sine data',c=col)
			ax = plt.gca()
			ax.set_xscale('log')
			
			if models is not None:
				if models == 'both' or models =='all':
					models = 'max hp'
				models = models.split()
				for model in models:
					self.SystemID(model)
					gain = 10**(self.dB[0]/20)*self.omegaN**2/np.sqrt((self.omegaN**2-omegas**2)**2+(2*self.zeta*self.omegaN*omegas)**2)
					plt.semilogx(omegas,20*np.log10(gain),label='model')

			plt.xlabel('Frequency [rad/s]')
			plt.ylabel('Amplitude [dB]')
			plt.legend()

		elif toPlot =='mag':
			plt.scatter(omegas/self.omegaN,dB-dB[0])
			ax = plt.gca()
			ax.set_xscale('log')
		elif toPlot == 'phase':
			plt.scatter(self.omegas,self.phase,label = 'Constant Sine data',c=col)
			ax = plt.gca()
			ax.set_xscale('log')	
			if models is not None:
				if models == 'both' or models == 'all':
					models = 'max hp'
				models = models.split()
				for model in models:
					self.SystemID(model)
					phase = (180/math.pi)*np.arctan2(-2*self.zeta*self.omegaN*omegas,self.omegaN**2-omegas**2)
					plt.semilogx(omegas,phase,label='model')
			plt.xlabel('Frequency [rad/s]')
			plt.ylabel('Phase [Rads]')
			plt.legend()

		elif toPlot == 'raw':
			data = np.genfromtxt(self.folder+'/'+self.datasets[0], delimiter=",",skip_header=1)
			t = data[:,0]
			pos = data[:,9]
			ce = data[:,4]
			plt.plot(t,pos, label="Displacement")
			plt.plot(t,ce,label = "Control Effort")
			plt.xlabel('Time [s]')
			plt.ylabel('Amplitude [counts]')
			plt.legend()

		
		ax = plt.gca()
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		if show:
			plt.show()
	def parameters(self,run):
		if run.mass == self.mass:
			raise ValueError("Mass of other run must be different!");
		if None in (run.omegaN,self.omegaN):
			raise ValueError("Must perform SystemID on both runs!")
		self.m = (run.mass*(run.omegaN**2)-self.mass*(self.omegaN**2))/((self.omegaN**2-run.omegaN**2))
		self.k = (self.m+self.mass)*self.omegaN**2;
		self.c = 2*self.zeta*np.sqrt((self.m+self.mass)*self.k)
		self.KHWOL = (10**(self.dB[0]/20))*self.k

	def print(self,method=None):
		print()
		print("Constant Sine: ", self.folder)
		if method is not None:
			label = self.SystemID(method = method, retLabel=True)
			print("ID method: ", label,"\n")
		if self.zeta is not None:
			print("omegaN = ", self.omegaN,"\n")
			print("zeta = ", self.zeta,"\n")
		if self.m is not None:
			print("m = ", self.m,"\n")
			print("k = ", self.k,"\n")
			print("c = ", self.c,"\n")
			print("KHWOL = ", self.KHWOL,"\n")

			
			

			