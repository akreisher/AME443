import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks
from scipy.optimize import curve_fit, fsolve
import os, sys, math

class StepResponse:
    def __init__(self,path,cutoff=None,tdiv=1):
        self.path = path;
        self.name = path.split('/')[-1]
        self.mass = ''
        self.string = ''
        self.stepInput = ''
        self.Kp = '0'
        self.Ki = '0'
        self.Kd = '0'

        argdict ={
            'm'  : 'mass',
            'r'  : 'stepInput',
            'Kp' : 'Kp',
            'KP' : 'Kp',
            'Ki' : 'Ki',
            'KI' : 'Ki',
            'Kd' : 'Kd',
            'KD' : 'Kd'
        }
        fileName = self.name
        #Parse experiment paramters from file name
        args = fileName.split("_")
        self.response = args[0]
        for arg in args[1:-1]:
            if arg[0] in argdict:
                setattr(self, argdict[arg[0]], arg[1:])
            if arg[:2] in argdict:
                setattr(self, argdict[arg[:2]], arg[2:])

        self.mass = float(self.mass)/1000
        self.springs = next(s for s in args if s[0] == 'k')
        self.stepInput = float(self.stepInput)

        self.Kp = float(self.Kp)
        self.Kd = float(self.Kd)
        self.Ki = float(self.Ki)

        self.runNum = int((args[-1].split("."))[0])
        #import all data into numpy array
        self.rawData = pd.read_csv(path)
        t = np.ndarray.flatten(self.rawData[["time"]].values)
        pos = np.ndarray.flatten(self.rawData[["actpos1"]].values)
        pks,_ = find_peaks(pos)
        #Check if data was already processed
        self.pPath = "/".join(path.split('/')[:-1])+"/PROCESSED/P"+fileName;

        if os.path.isfile(self.pPath):
            #Get start and end indices from file
            ind = np.genfromtxt(self.pPath,delimiter=",")
            end   = int(ind[1])
        elif cutoff is not None:
            end = cutoff;
        else:
            #Plot data and peaks
            plt.ion() #Interactive mode on
            plt.show()
            plt.plot(t,pos)
            plt.draw()
            plt.pause(0.001)

            #Prompt user for start and end peaks
            start = 1
            end = int(input("What time [s] to cutoff: "))-1

            #Close plot
            plt.ioff()
            Plt.close()

            #Store user input for future use
            with open(self.pPath,'w') as pfile:
                pfile.write(str(0)+","+str(end)+"\n")

        #Store cutoff data

        self.time = (t)[t<end]/tdiv

        self.pos = (pos)[0:len(self.time)]
        self.CE = np.ndarray.flatten(self.rawData[["iqcmd1"]].values)[0:len(self.time)]
        self.peaks = pks[pks<len(self.time)]

        #Default values
        self.omegaN = None
        self.zeta = None
        self.m = None
        self.K_HW = None

        pos = self.pos
        t = self.time
        pks = self.peaks

         #Steady State
        xInf = pos[-1]
        self.xInf = xInf
        self.ess =100* abs(xInf-self.stepInput)/self.stepInput

        # Overshoot
        xMax = max(pos)
        self.tPeak = t[np.argmax(pos)]
        self.overshoot = ((xMax-xInf)/xInf)*100;

        #rise time
        i = 0
        while pos[i]<xInf and i<len(t)-1:
            i=i+1
        #linear interpolation
        self.tRise = t[i-1]+(t[i]-t[i-1])*(xInf-pos[i-1])/(pos[i]-pos[i-1])

        #settling time
        i = len(t)-1
        while abs(pos[i]-xInf) < 0.02*xInf and i > 0:
            i = i-1
        self.tSettling = t[i]


    def SystemID(self, method='fit',toprint = False,retLabel = False):
        t = self.time
        pos = self.pos
        pks = self.peaks[:-2]
        if len(pks) > 3:
            pks = pks[:-1]
        xInf = self.xInf
        xMax = pos[pks[0]]
        label = None

        #Calculate omegaD from period
        T = np.mean(np.diff(t[pks[0:4]]))
        self.omegaD = 2*math.pi/(T)

        #Zeta by overshoot
        if (method=='overshoot'):

            #Calculate zeta
            logterm = np.log((xMax/xInf)-1)
            self.zeta = np.abs(logterm/np.sqrt(logterm**2+math.pi**2))

            self.omegaN = self.omegaD/math.sqrt(1-self.zeta**2)
            if toprint:
                print("System ID - Overshoot")
                print('x_max = ',xMax)
                print("omegaN = ",float(self.omegaN))
                print("zeta = ",float(self.zeta))
                print('\n')
            label = 'Overshoot'

        #Zeta by curvefittg
        elif (method == 'fit'):
            #Make mean(pos)=0
            pos = pos - xInf
            popt,pcov = curve_fit(lambda t,a,c:a*np.exp(c*t),t[pks],pos[pks])
            prod = -popt[1]
            #System to solve
            def syseqs(p):
                x ,y = p
                return [x*y-prod,x*np.sqrt(1-y**2)-self.omegaD]

            #Solve for params
            self.omegaN, self.zeta = fsolve(syseqs,(30,0.01))

            #Fix pos
            pos = pos+xInf

            if toprint:
                print("System ID - Fit")
                print("zeta = ",float(self.zeta))
                print("omegaN = ",float(self.omegaN))
                print('\n')

            label = "Curve Fit"


        #Zeta by rise time
        elif (method=='rise'):
            #find rise time



            phi = math.pi-tRise*self.omegaD
            self.zeta = fsolve(lambda x:np.arctan(np.sqrt(1-x**2)/x)-phi,0.01)
            self.omegaN = self.omegaD/np.sqrt(1-self.zeta**2)

            if toprint:
                print("System ID - Rise")
                print("t_rise = ",tRise)
                print("zeta = ",float(self.zeta))
                print("omegaN = ",float(self.omegaN))
                print('\n')
            label = 'Rise time'

        else:
            print("Method must be overshoot or rise")
        if retLabel:
            return label

    #Plot, can supply args in dict
    def plot(self,fig=None, data = True, models = None, show=True, peaks = False):
        #Get figure and axis
        fig = plt.figure(fig)
        fig.patch.set_facecolor('white')
        ax = plt.gca()
        t=self.time
        pos = self.pos
        pks = self.peaks
        #Plot data (default if no args given)
        if data:
            plt.plot(t,pos,label="Experimental data")
            if peaks:
                plt.scatter(t[pks],pos[pks])

        #Plot model
        if models is not None:
            #Plot one model
            if models == 'all':
                models = 'overshoot rise fit'
            models = models.split(" ")
            for model in models:
                label = self.SystemID(method=model,retLabel = True)
                if label is not None:
                    phi = math.atan(np.sqrt(1-self.zeta**2)/self.zeta)
                    modelCurve = self.xInf*(1-(np.exp(-self.zeta*self.omegaN*t)*np.sin(self.omegaD*t+phi)/np.sqrt(1-self.zeta**2)))
                    plt.plot(t,modelCurve,label = label+" model")

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


    def parameters(self,run):
        if run.mass == self.mass:
            raise ValueError("Mass of other run must be different!");
        if None in (run.omegaN,self.omegaN):
            raise ValueError("Must perform SystemID on both runs!")
        self.m = (run.mass*(run.omegaN**2)-self.mass*(self.omegaN**2))/((self.omegaN**2-run.omegaN**2))
        self.k = (self.m+self.mass)*(self.omegaN**2)
        self.c = 2*self.zeta*np.sqrt((self.m+self.mass)*self.k)
        self.K_HW = self.k*self.xInf/self.stepInput

    #Clear processed data for reprocessing
    def clear(self):
        os.remove(self.pPath);

    def print(self,method=None):
        print()
        print("Step Response: ", self.file)
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
            print("KHW = ", self.K_HW,"\n")
