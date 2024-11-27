import numpy as np
from scipy.integrate import ode
import matplotlib.pylab as plt
import matplotlib
from matplotlib import ticker
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as colors
from matplotlib.legend_handler import HandlerTuple
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 22}) # update fontsize

class particle():
    def __init__(self, ChargedParticleFlag, x0Vec, v0Vec, E, B):

        # Assign constants depending on charged particle chosen
        if ChargedParticleFlag == 0:
            self.q = 1.6022e-19
            self.m = 1.6276e-27
        elif ChargedParticleFlag == 1:
            self.q = -1.6022e-19
            self.m = 9.1094e-31
        else:
            assert ChargedParticleFlag == 0 or ChargedParticleFlag == 1, "particle identifier flag chosen incorrectly. Should be 0 for ions and 1 for electrons!"
        self.x0Vec = x0Vec
        self.v0Vec = v0Vec

        self.EfieldVec = E
        self.BfieldVec = B

    def Get_dvdt(self,vVec):
        return self.q/self.m*(self.EfieldVec + np.cross(vVec,self.BfieldVec))

    def Get_dxdt(self,vVec):
        return vVec


    def rk4(self,t,dt,xVec,vVec):
        
        k1 = dt * self.Get_dxdt(vVec)
        h1 = dt * self.Get_dvdt(vVec)
        
        k2 = dt * self.Get_dxdt(vVec + h1/2)
        h2 = dt * self.Get_dvdt(vVec + h1/2)

        k3 = dt * self.Get_dxdt(vVec + h2/2)
        h3 = dt * self.Get_dvdt(vVec + h2/2)

        k4 = dt * self.Get_dxdt(vVec + h3)
        h4 = dt * self.Get_dvdt(vVec + h3)

        xVec = xVec + 1/6 * (k1 + 2*k2 + 2*k3 + k4)
        vVec = vVec + 1/6 * (h1 + 2*h2 + 2*h3 + h4)
        t = t + dt

        return t, xVec, vVec


    def pushParticle(self, dt, numsteps):
        
        t0 = 0
        count = 0
        xVec = self.x0Vec
        vVec = self.v0Vec
        tVec = t0

        tVecSave = [tVec]
        vVecSave = [vVec]
        xVecSave = [xVec]
        
        while count < numsteps:
            t,x,v = self.rk4(tVec, dt, xVec, vVec)
            
            tVec = t
            xVec = x
            vVec = v

            tVecSave.append(tVec)
            xVecSave.append(xVec)
            vVecSave.append(vVec)

            count += 1

        vVecData = np.zeros((3,len(tVecSave)))
        xVecData = np.zeros((3,len(tVecSave)))

        for i in range(len(tVecSave)):
            vVecData[:,i] = vVecSave[i][:]
            xVecData[:,i] = xVecSave[i][:]

        return tVecSave, xVecData, vVecData

    def PlotParticle(self, t, x, v):
        
        fig = plt.figure(clear=True)
        fig.set_tight_layout(True)
        ax = fig.add_subplot(111,projection='3d')

        ax.plot3D(x[:,0], x[:,1], x[:,2],
                 '-k', linewidth=2)
        plt.show()
        

            
        


