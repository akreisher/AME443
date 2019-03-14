import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit, fsolve
import os, sys, math

"""
AUTHOR: Alex Kreisher
AME 443

Class representing a data run for the Free Response of a 
second order underdamped system. Perfroms System ID to 
calculate characteristic params. based on log. decrement or
curve fitting methods.

"""

class FreeResponse:
	def __init__(self,fileName,weekNum):
		self.weekNum = weekNum;
		self.file = fileName;
		self.path = "Week "+str(weekNum)+"/"+fileName;
		#Parse experiment paramters from file name
		args = fileName.split("_");
		self.mass = float((args[1])[1:len(args[1])])/1000;
		self.springs = args[2];
		self.runNum = int((args[3].split("."))[0]);
		
		#import all data into numpy array
		self.rawData = np.genfromtxt(self.path, delimiter=",",skip_header=1)
		self.peaks,_ = find_peaks(self.rawData[:,9]);
		#Check if data was already processed
		self.pPath = "Week "+str(weekNum)+"/PROCESSED/P"+fileName;
		if os.path.isfile(self.pPath):
			#Get start and end indices from file
			ind = np.genfromtxt(self.pPath,delimiter=",");
			self.start = int(ind[0]);
			self.end   = int(ind[1]);
		else:
			#Plot data and peaks
			plt.ion(); #Interactive mode on
			plt.show();
			plt.title(self.file)
			plt.plot(self.rawData[:,0],self.rawData[:,9])
			plt.scatter(self.rawData[:,0][self.peaks],self.rawData[:,9][self.peaks],c='r')
			plt.draw();
			plt.pause(0.001)

			#Prompt user for start and end peaks
			self.start = int(input("What peak to start at? (Starting at 1): "));
			self.end = int(input("How many peaks to ignore? (Starting at 1 from end): "));
			
			#Close plot
			plt.ioff();
			plt.close()

			#Store user input for future use
			with open(self.pPath,'w') as pfile:
				pfile.write(str(self.start)+","+str(self.end)+"\n")
		
		#Store cutoff data
		self.time = (self.rawData[:,0])[(self.peaks[self.start-1]):]/10
		self.time = self.time - self.time[0];
		self.pos = (self.rawData[:,9])[self.peaks[self.start-1]:len(self.rawData[:,9])]
		self.peaks = self.peaks[self.start-1:-self.end]-self.peaks[self.start-1];
		
		#Default values
		self.omegaN = None;
		self.zeta = None;
		self.m = None;
		self.c = None;
		self.k = None;

	def SystemID(self, method='fit',toPrint = False,retLabel = False):
		self.x0 = self.pos[0];
		label = None
		#Calculate omegaD from period
		T = np.empty([len(self.peaks)-1,1]);
		for n,ind in enumerate(self.peaks[1:len(self.peaks)]):
				T[n] = self.time[ind]/(n+1)
		T = np.mean(T);
		self.omegaD = 2*math.pi/(T);

		#Zeta by logarithmic decrement
		if (method=='log'):
			#Set up empty arrays
			zetas = np.empty([len(self.peaks)-1,1]);
			#Calculate decrement for each peak
			for n,ind in enumerate(self.peaks[1:len(self.peaks)]):
				xn = self.pos[ind]
				zetas[n] = math.log(self.x0/xn)/math.sqrt(4*(math.pi**2)*((n+1)**2)+math.log(self.x0/xn)**2);
			#Take average zeta
			self.zeta = np.mean(zetas);
			self.omegaN = self.omegaD/math.sqrt(1-self.zeta**2);
			label = 'Logarithmic decrement'
		#Zeta by curvefitting
		elif (method=='fit'):
			#Fit peaks to a*e^(bt)
			popt,pcov = curve_fit(lambda t,a,b:a*np.exp(b*t),self.time[self.peaks],self.pos[self.peaks]) 
			prod = -popt[1]; #b
			#System of equations to be solved
			def syseqs(p):
				x,y = p
				return [x*y-prod,x*np.sqrt(1-y**2)-self.omegaD]
			#Solve for omegaN and zeta
			self.omegaN, self.zeta = fsolve(syseqs,(30,0.01))
			label = 'Curve fit'
		else:
			#Bad method
			print("System ID method must be log or fit")
		if retLabel:
			return label

	#Plot, can supply args in dict or as keywords
	#Args: 
		#fig:figure number
		#data: show processed data
		#models: mdoels to show as string w mehtods as words seperated by space
		#show: show plot
		#peaks: show peaks
	def plot(self,fig=None,data = True, models = None,show=True,peaks = False):
		#Get figure and axis
		fig = plt.figure(fig)
		ax = plt.gca()
		t = self.time
		pos = self.pos
		pks = self.peaks
		#Plot data (default if no args given)
		if data:
			plt.plot(t,pos,label="Experimental data")
			if peaks:
				plt.scatter(t[pks],pos[pks])
		
		#Plot model
		if models is not None:
			if models == 'all':
				models = 'log fit'
			models = models.split(" ")
			for model in models:
				label = self.SystemID(method=model, retLabel=True)
				if label is not None:
					beta = math.atan(self.zeta/math.sqrt(1-self.zeta**2))
					modelCurve = (self.x0/np.sqrt(1-self.zeta**2))*np.cos(self.omegaD*self.time-beta)*np.exp(-self.zeta*self.omegaN*self.time)
					plt.plot(self.time,modelCurve,label = label+' model');
				
		#Format plot
		plt.xlim(0,self.time[len(self.time)-1])
		plt.xlabel('Time [s]')
		plt.ylabel('Displacement [counts]')
		plt.legend()
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		#Show
		if show:
			plt.show()
	
	#Calculate physical paramters based on 2 runs w/ diff. masses
	def parameters(self,run):
		#Error if same mass
		if run.mass == self.mass:
			raise ValueError("Mass of other run must be different!");
		#No values for characteristic params
		if None in (run.omegaN,self.omegaN):
			raise ValueError("Must perform SystemID on both runs!")
		
		#Calc m, k, c
		self.m = (run.mass*(run.omegaN**2)-self.mass*(self.omegaN**2))/((self.omegaN**2-run.omegaN**2))
		self.k = (self.m+self.mass)*self.omegaN**2;
		self.c = 2*self.zeta*np.sqrt((self.m+self.mass)*self.k)
	
	#Clear processed data for reprocessing
	def clear(self):
		os.remove(self.pPath);

	def print(self,method=None):
		print()
		print("Free Response: ", self.file)
		if method is not None:
			label = self.SystemID(method = method, retLabel=True)
			print("ID method: ", label, "\n")
		if self.zeta is not None:
			print("omegaN = ", self.omegaN,"\n")
			print("zeta = ", self.zeta,"\n")
		if self.m is not None:
			print("m = ", self.m,"\n")
			print("k = ", self.k,"\n")
			print("c = ", self.c,"\n")

