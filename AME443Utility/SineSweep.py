import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import os, sys, math

class SineSweep:
	def __init__(self,fileName,weekNum):
		self.weekNum = weekNum;
		self.file = fileName;
		self.path = "Week "+str(weekNum)+"/"+fileName;
		#Parse experiment paramters from file name
		args = fileName.split("_");
		self.mass = float((args[1])[1:len(args[1])])/1000;
		self.springs = args[2];
		self.inputVoltage = float((args[3])[1:len(args[3])])
		self.runNum = int((args[5].split("."))[0])
		
		#import all data into numpy array
		self.rawData = np.genfromtxt(self.path, delimiter=",",skip_header=1)
		self.time = self.rawData[:,0]
		self.pos = self.rawData[:,9]
		self.ce = self.rawData[:,4]
		t = self.time
		pos = self.pos
		ce = self.ce
		pks,_ = find_peaks(ce)
		Ts = np.diff(t[pks])
		omegas = 2*math.pi/Ts
		self.rawOmegas = omegas;

		popt,pcov = curve_fit(lambda t,a,b:a*t+b,t[pks[1:len(pks)]],omegas)
		self.popt = popt

		pks,_ = find_peaks(pos)
		self.dB = 20*np.log10(np.abs(pos[pks]/self.inputVoltage))
		self.omegas = popt[0]*t[pks]+popt[1]
		self.zeta = None
		self.m = None




	def SystemID(self,method='max',retLabel=False):
		
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
			label='Max magnification'
		else:
			print('System ID method must be hp or max')
		if retLabel:
			return label



	def plot(self,fig=None,show=True, toPlot='gain',models=None):
		omegas = self.omegas
		dB=self.dB
		fig = plt.figure(fig)
		if toPlot=='gain':
			plt.semilogx(omegas,dB,label='Sine Sweep data')
			
			if models is not None:
				if models == 'both' or models == 'all':
					models = 'max hp'
				models = models.split(" ")
				for model in models:
					label = self.SystemID(method=model,retLabel=True)
					if label is not None:
						gain = 10**(self.dB[0]/20)*self.omegaN**2/np.sqrt((self.omegaN**2-omegas**2)**2+(2*self.zeta*self.omegaN*omegas)**2)
						plt.semilogx(omegas,20*np.log10(gain),label=label+' model')
			
			plt.xlabel('Frequency [rad/s]')
			plt.ylabel('Gain [dB]')
		if toPlot =='mag':
			plt.semilogx(omegas/self.omegaN,dB-dB[0])
			plt.xlabel(r'$\omega/{\omega}_n$ [-]')
			plt.ylabel(r'|M(i$\omega$)| [dB]')
		if toPlot == 'raw':
			f, (ax1,ax2) = plt.subplots(2,1,sharex=True,sharey=True)
			ax1.plot(self.time,self.pos, label= "Displacement")
			ax1.spines['top'].set_visible(False)
			ax1.spines['right'].set_visible(False)
			ax1.spines['bottom'].set_visible(False)
			ax1.legend()
			ax1.tick_params(axis='x',which='both', bottom=False,labelbottom=False)
			ax2.plot(self.time,self.ce, label = "Control Effort", color = 'orange')
			plt.xlabel('Time [s]')
			f.text(0.02, 0.5, 'Amplitude [counts]', ha='center', va='center', rotation='vertical')

		if toPlot == 'calib':
			pks,_ = find_peaks(self.ce)
			plt.plot(self.time[pks[1:]],self.rawOmegas,label = "Raw data")
			pks,_ = find_peaks(self.pos)
			plt.plot(self.time[pks],self.omegas, label="Calibration Curve")

			plt.xlabel('Time [s]')
			plt.ylabel('Frequency [rad/s]')


		
		plt.legend()
		ax = plt.gca()
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		if show:
			plt.show()
			plt.savefig('plot')
			plt.clf()
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
		print("Sine Sweep: ", self.file)
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
