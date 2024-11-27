import ChargedParticlePusher as CPP
import numpy as np

ChargedParticleFlag = 1
x0Vec = np.array([0,0,0])
v0Vec = np.array([1,1,0])
E = np.array([0,0,0])
B = np.array([0,0,5.3])

MyParticle = CPP.particle(ChargedParticleFlag,
                         x0Vec, v0Vec,
                         E, B)
dt = 1.e-8
numsteps = 50
t, x, v = MyParticle.pushParticle(dt,numsteps)
MyParticle.PlotParticle(t,x,v)
