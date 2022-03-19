import matplotlib.pyplot as plt
#from numpy import *
import numpy as np
import os
import glob
from multiprocessing import Pool
from CIFTVelocityLargeRB import loadCIFTFile
import lmfit
from FittingTools import *





class CIFTBFieldLargeRB:

  def __init__(self, directoryNamebField, directoryNameGeometry,
               nRings, nSensorsAtRing, shiftRingSensors = None):
    """ Loads the magnetic field data of the experiment and 
        loads the sensor geometrical information. 

        ATTENTION: The program assumes that the sensors along one 
                   rail over height is saved consecutively.

       Arguments
       directoryNamebField   -- path to the bfield data
       directoryNameGeometry -- path to the geometry info
       nRings                -- number of sensor rings over height
       nSensorsAtRing        -- number of sensors at one ring
       shiftRingSensors -- Integer-array. The number of elements
                           is the same as the number of rings.
                           It defines for each ring the number of 
                           shifting of the sensors so that the 
                           first sensor has the smalles position angle.
                           Needed for sinus fit
    """
    self.nRings         = nRings
    self.nSensorsAtRing = nSensorsAtRing 
    bfieldFileNames =  np.sort(glob.glob(directoryNamebField + "/*.dat"))
    
    self.times = np.array(list(map(lambda s:int(os.path.basename(s).
                                              replace(".dat","").
                                                split("_")[-1]),
                                   bfieldFileNames)))
    p = Pool()
    self.bList    = list(p.map(loadCIFTFile, bfieldFileNames))
    self.b        = np.row_stack(self.bList)
    
    self.pSens = np.loadtxt(directoryNameGeometry+"/sensorPos.dat", skiprows=1)
    self.nSens = np.loadtxt(directoryNameGeometry+"/sensorDir.dat", skiprows=1)

    self.shiftRingSensors = shiftRingSensors
      

    
  def plotSpatialTemporalVertical(self, rangeInnT, nCols = None,
                                  figNr = 2):
    """ Generates a plot for each sensor column visualising the 
        magnetic field magnitude over time and height for each sensor rail. 
        Arguments
        rangeInnT -- list of the form [minRange, maxRange] in nT
        nCols     -- consider only the first nCols
        figNr     -- number of the figure
    """
    
    # 1.) Get sensor grid
    z = self.pSens[0:self.nRings,2]
    X,Y = np.meshgrid(self.times, z)

    # 2.) get levels
    lev = np.linspace(rangeInnT[0], rangeInnT[1], (rangeInnT[1] - rangeInnT[0]))

    # 3.) get number of plots:
    if nCols is None:
      N = int(np.ceil(self.pSens.shape[0] / self.nRings))
    else:
      N = nCols
      
    # 4.) Plot values
    fig, axes = plt.subplots(nrows=N, sharex=True, num=figNr, clear=True)
    fig.set_size_inches([16, max(8,2.2*N)])
    for i in range(N):
      data = np.transpose(self.b[:,
                                 self.nRings*i:
                                 min(self.nRings*(i+1),
                                     self.b.shape[1])])*1e9
      pl = axes[i].contourf(X,Y, data,
                            cmap=plt.cm.jet,
                            levels = lev)
      
      axes[i].set_ylabel("col %i\n z in mm" % (i+1))

      # Colorbar
      cb = fig.colorbar(pl, ax = axes[i])
      cb.set_ticks([lev.min(), 0, lev.max()])

    axes[-1].set_xlabel("time in s")
    fig.tight_layout()



  def plotSpatialTemporalHorizontal(self, rangeInnT,
                                    figNr = 3):
    """ Generates a plot for each sensor ring visualising the 
        magnetic field magnitude over time and azimiuth for each sensor ring. 
        Arguments
        rangeInnT -- list of the form [minRange, maxRange] in nT
        figNr     -- number of the figure
    """
    if self.shiftRingSensors is None:
      self.shiftRingSensors = np.zeros(self.nRings, dtype=int)

    # 1.) Get sensor grid
    azimut = np.arctan2(self.pSens[:,0], self.pSens[:,1]) *180/np.pi
    z = self.pSens[0:self.nRings,2]

    # 1a.) List of indices per ring
    L = [[] for i in range(z.shape[0])]
    for i in range(self.pSens.shape[0]):
      L[i % self.nRings] += [i]

    
    # 2.) get levels
    lev = np.linspace(rangeInnT[0], rangeInnT[1], (rangeInnT[1] - rangeInnT[0]))

    # 3.) get number of plots:
    N = len(L)
      
    # 4.) Plot values
    fig, axes = plt.subplots(nrows=N, sharex=True, num=figNr, clear=True)
    fig.set_size_inches([16, max(8,2.2*N)])
    for i in range(N):
      X,Y = np.meshgrid(self.times, azimut[L[i]])
      idxList = np.roll(L[i],    self.shiftRingSensors[i])
      pl = axes[i].contourf(X,Y, np.transpose(self.b[:,idxList]) * 1e9,
                            cmap = plt.cm.jet,
                            levels = lev)
      
      axes[i].set_ylabel("ring %i\nangle [°]" % (i+1))
      axes[i].set_ylim(-180,180)
      
      # Colorbar
      cb = fig.colorbar(pl, ax = axes[i])
      cb.set_ticks([lev.min(), 0, lev.max()])

    axes[-1].set_xlabel("time in s")
    fig.tight_layout()



    
  def plotSensorsOverHeight(self, t, rangeInnT, nCols = None):
    """ Plots the magnetic field data over the height for the given time 
        instance """
    # 1.) Plot title
    plt.clf()
    plt.title("Time %d s " % (t))
    
    # 2.) get number of plots:
    if nCols is None:
      N = int(np.ceil(self.pSens.shape[0] / self.nRings))
    else:
      N = nCols
    Nmax = self.pSens.shape[0]
      
    # 3.) Plot sensors
    for i in range(N):
      plt.plot(self.pSens[self.nRings*i:min(Nmax,self.nRings*(i+1)), 2],
               self.b    [t, self.nRings*i:min(Nmax,self.nRings*(i+1))] * 1e9,
               'o-', label="%d-%d" % (N*i,min(Nmax,N*(i+1)))) 
    plt.legend()
    plt.grid()
    
    # 4.) Annotate sensors
    for i in range(self.pSens.shape[0]):
      plt.annotate("%d" % i, (self.pSens[i, 2],  self.b[t, i] * 1e9))

    
    

  def plotSinusFitAtRings(self, t, rangeInnT, figNr = 4):
    """ Plots the measured magnetic field for one time instance 
        in dependence of its angular position. 
    """
    if self.shiftRingSensors is None:
      self.shiftRingSensors = np.zeros(self.nRings, dtype=int)
      
    # 1.) Get sensor grid
    azimut = np.arctan2(self.pSens[:,0], self.pSens[:,1]) *180/np.pi
    z = self.pSens[0:self.nRings,2]

    # 1a.) List of indices per ring
    L = [[] for i in range(z.shape[0])]
    for i in range(self.pSens.shape[0]):
      L[i % self.nRings] += [i]

    N = len(L)
      
    # 2.) Init figure
    fig, axes = plt.subplots(nrows=N, sharex=True, num = figNr, clear=True)
    fig.set_size_inches([16, 11])

    # Set angles and color cycle
    allAngles = np.linspace(-180, 180, 181)
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # 3.) Generate the plot
    for i in range(N):
        x = np.roll(azimut[L[i]],    self.shiftRingSensors[i])
        y = np.roll(self.b[t, L[i]], self.shiftRingSensors[i])
        p = SinusFitDegree(x,y)
        c = cycle[i % len(cycle)]
        axes[i].plot(x, y * 1e9, "--o", c=c, label="ring %d" % (i+1))
        axes[i].plot(allAngles, p.calc(allAngles) * 1e9, "-", c=c)
        axes[i].set_ylabel("ring %d\nb [nT]" % (i+1))
        axes[i].set_ylim(rangeInnT[0], rangeInnT[1])  
        
    axes[-1].set_xlabel("angle in deg")


    
  def getAmplAndAngleAtRings(self):
    """ Returns the angle and the amplitude of the sinus fit for 
        each sensor over time """

    if self.shiftRingSensors is None:
      self.shiftRingSensors = np.zeros(self.nRings, dtype=int)
      
    # 1.) Get sensor grid
    azimut = np.arctan2(self.pSens[:,0], self.pSens[:,1]) *180/np.pi
    z = self.pSens[0:self.nRings,2]

    # 1a.) List of indices per ring
    L = [[] for i in range(z.shape[0])]
    for i in range(self.pSens.shape[0]):
      L[i % self.nRings] += [i]

    N = len(L)    

    self.ringAngle = np.zeros((self.times.shape[0], N))
    self.ringAmp   = np.zeros_like(self.ringAngle)
    for tid in range(self.times.shape[0]):
      for i in range(N):
        x = np.roll(azimut[L[i]],    self.shiftRingSensors[i])
        y = np.roll(self.b[tid, L[i]], self.shiftRingSensors[i])
        p = SinusFitDegree(x,y)
        
        self.ringAngle[tid, i] = p.getFitParameter()[1]
        self.ringAmp[tid, i]   = p.getFitParameter()[0]

    return self.ringAngle, self.ringAmp
