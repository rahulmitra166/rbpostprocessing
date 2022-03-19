import numpy as np
import matplotlib.pyplot as plt
from numpy import *
from matplotlib.pyplot import *
import glob
import os
from multiprocessing import Pool 
from functools import reduce
        

def loadCIFTFile(s):
    """ Checks if a file exists and reads it ASCII file """
    if os.path.isfile(s) and os.access(s, os.R_OK):
        return np.loadtxt(s, skiprows=1)
    else:
        print("File %s does not exists or is not readable" % s)
        return None

    
def loadListOfFiles(directoryName, fileNamePattern,
                    dontUseBinaryFiles = False):
    """ Loads a list of files in memory. The time is given in the 
        file name  in seconds with the format <prefix>_<time>.<ext>. 
        From this file list a binary representation of the data 
        is given in the directory:
        data.npy -> the contents of the files
        dataTime -> the time stamps.
    
        Later the binary files can be directly read speeding up the 
        process.
    """

    compressedFiles = "%s/data.npy" % directoryName
    timeFiles       = "%s/dataTime.npy" % directoryName

    if os.path.isfile(compressedFiles) and os.access(compressedFiles, os.R_OK)\
       and os.path.isfile(timeFiles) and  os.access(timeFiles, os.R_OK) \
       and not dontUseBinaryFiles:

        print("Load compacted binary data")

        # If compressed file exists
        D     = np.load(compressedFiles)
        times = np.load(timeFiles)

    else:
        print("Load file list")
        
        # 1.) Get file names
        fileNames = np.sort(glob.glob("%s/%s" % (directoryName,
                                                 fileNamePattern)))
    
        # 2.) Parallel read
        p = Pool()
        D = np.array(list(p.map(loadCIFTFile, fileNames)))
    
    
        # 3.) Generate times from file names
        times = np.array(list(map(lambda s:int(os.path.basename(s).
                                               replace(".dat","").
                                               split("_")[-1]),
                                  fileNames)))

        # Save binary file
        if not  dontUseBinaryFiles:
            np.save(compressedFiles, D)
            np.save(timeFiles, times)

    return times, D

    
    
class CIFTVelocityLargeRB:
    
    def __init__(self, directoryNameVelo, volumePointsDir):
        
        # 1.) Load velocity data
        self.time, self.velo = loadListOfFiles(directoryNameVelo + "/",
                                                "*.dat")
        
        # 2.) Load position data & calculate the radius
        self.pVol  = np.loadtxt(volumePointsDir + "/points_vol.dat", skiprows=1)
        self.pR    = np.sqrt(self.pVol[:,0]**2 + self.pVol[:,1]**2)
        
        # 3.) Reynolds number factor for velocity
        self.nu      = 3.44e-7 # m^2/s (GaInSn)
        self.charLen = 0.64
        self.ReFac   = self.charLen / self.nu
        
        self.__initSpatialInfo__()
        self.__calc__()
        

        
    def __initSpatialInfo__(self):
        """ Extracts the z-Levels for slicing"""
        self.zLevels = np.sort(list(set(self.pVol[:,2].round(4))))
        self.zLevelIndices = [[] for i in range(self.zLevels.shape[0])]
        for i in range(self.pVol.shape[0]):
            diff = abs(self.pVol[i,2] - self.zLevels)
            j    = argmin(diff)
            self.zLevelIndices[j] += [i] 


            
    def __calc__(self):
        """ Calculates important parameters for every time step:
             
            1.) mean, max, std of the velocity fielf
            2.) Energy of vz component
            3.) Direction, angle and magnitude of the 
                horizontal velocity at each z-level
        """ 
        # Velocity magnitude
        self.Magnitude = np.sqrt(self.velo[:,:,0]**2 + \
                                 self.velo[:,:,1]**2 + \
                                 self.velo[:,:,2]**2)

        
        # Mean and max velocity magnitude
        self.MagMean = self.Magnitude.mean(axis=1)
        self.MagMax  = self.Magnitude.max(axis=1)
        self.MagStd  = self.Magnitude.std(axis=1)

        # Energy vz
        self.sumVz = np.zeros((self.time.shape[0],
                               self.zLevels.shape[0]))
        # Absolute mean vz
        self.absMeanVz = np.zeros((self.time.shape[0],
                                self.zLevels.shape[0]))

        # Mean velocity all
        self.meanVall = np.zeros((self.time.shape[0],
                                  self.zLevels.shape[0]))

        # Extract the data
        for k in range(self.time.shape[0]):
            for i in range(self.zLevels.shape[0]):

                # sum(vz^2) per level  
                res = (self.velo[k, self.zLevelIndices[i], 2]**2).sum()
                self.sumVz[k,i] = res

                # avg(abs(vz)) per level  
                res = np.abs(self.velo[k, self.zLevelIndices[i], 2]).mean()
                self.absMeanVz[k,i] = res

                # get all velocities per level
                res = self.Magnitude[k,  self.zLevelIndices[i]].mean()
                self.meanVall[k,i] = res
                
        #self.extractDirectionVectorAtZLevels()

        
    def extractDirectionVectorAtZLevels(self, r = 1):
        # Direction vx,vy
        self.dirVxVy = np.zeros((self.time.shape[0],
                                 self.zLevels.shape[0], 2))

        for k in range(self.time.shape[0]):
            for i in range(self.zLevels.shape[0]):
                v = self.velo[k, self.zLevelIndices[i], 0:2]
                idxL = np.where(self.pR[self.zLevelIndices[i]] < r)
                self.dirVxVy[k, i, :] = v[idxL[0], :].mean(axis=0)
                
        self.angleVxVyRaw = np.zeros((self.time.shape[0],self.zLevels.shape[0]))
        self.angleVxVy    = np.zeros((self.time.shape[0],self.zLevels.shape[0]))
        self.magVxVy      = np.zeros((self.time.shape[0],self.zLevels.shape[0]))
                       
        for k in range(self.time.shape[0]):
            for i in range(self.zLevels.shape[0]):
                self.angleVxVy[k,:] = np.unwrap(np.arctan2(self.dirVxVy[k,:,1],
                                                           self.dirVxVy[k,:,0]))
                self.angleVxVyRaw[k,:] = np.arctan2(self.dirVxVy[k,:,1],
                                                    self.dirVxVy[k,:,0])
                self.magVxVy[k,:] = np.sqrt(self.dirVxVy[k,:,0]**2 +
                                            self.dirVxVy[k,:,1]**2)


                

    #####################################################################
    #
    # PLOTTING
    #
    ####################################################################
    def plotSpatiallyAveragedVelocity(self, asReynoldsNum):
        """ Plots the spatially averaged velocity 
           over time in m/s. If asReynoldsNumber = True, 
           as Reynolds nummer"""
        if asReynoldsNum:
            plot(self.time, self.MagMean*self.ReFac, label="mean")
            ylabel("Re")
        else:
            plot(self.time, self.MagMean, label="mean")
            ylabel("v in m/s")
        xlabel("t in s")


        
    def plotSpatiallyMaxVelocity(self, asReynoldsNum):
        """ Plots the spatially maximum velocity 
           over time in m/s. If asReynoldsNumber = True, 
           as Reynolds nummer"""
        if asReynoldsNum:
            plot(self.time, self.MagMax*self.ReFac, label="max")
            ylabel("Re")
        else:
            plot(self.time, self.MagMax, label="max")
            ylabel("v in m/s")
        xlabel("t in s")


        
    def getSpatiallyMeanMaxVelocityImageSeq(self, prefixFileName,
                                            startTime, endTime,
                                            deltaT,
                                            asReynoldsNum = True,
                                            plotMaxVelocity = False):
        """ Generate a list of images with the given prefixFileName 
          for a time intervall given by startTime and endTime in sec
          with the method plotSpatiallyAveragedVelocity and 
          plotSpatiallyMaxVelocity. The actual time is given by
          by a vertical bar and the interval is actual time +- deltaT.
        """
        clf()

        # Plot averaged velocity
        self.plotSpatiallyAveragedVelocity(asReynoldsNum)

        # Plot max velocity
        if plotMaxVelocity:
            self.plotSpatiallyMaxVelocity(asReynoldsNum)

        # Marker of the current position
        verticalLine = axvline(startTime)

        for k in range(startTime, endTime):
            title("time: %d s" % k)
            verticalLine.set_data([k,k], [0,1])
            xlim(k - deltaT, k + deltaT)
            draw()
            savefig("%s_%06d.png" % (prefixFileName, k))
        
        
        
    def saveSpatiallyAveragedData(self, fileName):
        A = np.column_stack((self.time,
                             self.MagMean, self.MagStd, self.MagMax,
                             self.MagMean*self.ReFac, self.MagStd*self.ReFac,
                             self.MagMax*self.ReFac))
        headerStr = "time in s, v mean in m/s, std in m/s, max in m/s, " + \
                    "Re mean, Re std, Re max"
        savetxt(fileName, A, header=headerStr)

        
        
    def plotVectorVxVyAtLevel(self, timeID, levelID):
        """ Plots the vx vy vectors at a given levelID"""
        quiver(self.pVol[self.zLevelIndices[levelID],0],
               self.pVol[self.zLevelIndices[levelID],1],
               self.velo[timeID,self.zLevelIndices[levelID],0],
               self.velo[timeID,self.zLevelIndices[levelID],1])      
        xlabel("x in m")
        ylabel("y in m")
        gca().set_aspect("equal")

        
    def plotVectorVxVyTopBottom(self,timeID):
        """ Plots the vx vy vectors at top and bottom"""
        clf()
        subplot(2,1,1)
        title("bottom")
        self.plotVectorVxVyAtLevel(timeID,1)
        xlabel("")
        for k in gca().get_xticklabels():
            k.set_visible(False)
        
        subplot(2,1,2)
        title("top")
        N = len(self.zLevels)
        self.plotVectorVxVyAtLevel(timeID,N-2)
    


    def getVectorVxVyTopBottomImageSeq(self, prefixFileName,
                                       startTime, endTime):
        """ Generate a list of images with the given prefixFileName 
          for a time intervall given by startTime and endTime in sec
          with the method plotVectorVxVyTopBottom."""
        
        for k in range(startTime, endTime):
            clf()
            title("time: %d s" % k)
            idx = searchsorted(self.time, k)
            self.plotVectorVxVyTopBottom(idx) 
            savefig("%s_%06d.png" % (prefixFileName, k))
                        
                         
        
    def plotAngleTopBottomOverTime(self):
        """ Plots the angle of the spatial average of the horizontal velocity
           at top and bottom """
        N = len(self.zLevels)
        angleBot = np.rad2deg(self.angleVxVy[:,1])
        angleTop = np.rad2deg(self.angleVxVy[:,N-2])
        plot(self.time, angleBot, label="bot")
        plot(self.time, angleTop, label="top")
        legend()
        xlabel("t in s")
        ylabel("angle in degree")
        

        
    def plotAngleAtLevelOverTime(self, levelID):
        """ Plots the angle of the spatial average of the horizontal velocity
           at the given levelID """
        angle = np.rad2deg(self.angleVxVy[:, levelID])
        p = plot(self.time, angle,
                 label="level %0.3f m" % self.zLevels[levelID])
        xlabel("t in s")
        ylabel("angle in degree")
        return p


    
    def plotMagTopBottomOverTime(self):
        """ Plots the magnitude of the spatial average of the horizontal 
        velocity at top and bottom """
        N = len(self.zLevels)
        magBot = self.magVxVy[:,1]
        magTop = self.magVxVy[:,N-2]
        plot(self.time, magBot, label="bot")
        plot(self.time, magTop, label="top")
        legend()
        xlabel("t in s")
        ylabel("v in m/s")


        
    def plotMagAtLevelOverTime(self, levelID):
        """ Plots the magnitude of the spatial average of the horizontal 
        velocity at the given levelID """
        mag = self.magVxVy[:, levelID]
        p=plot(self.time, mag, label="level %0.3f m" % self.zLevels[levelID])
        xlabel("t in s")
        ylabel("v in m/s")
        return p
       

    def plotAngleMagOverHeight(self, idxTime, plotAngle = True,
                               angleRange = None, plotMag = True,
                               magVeloRange = None):
        """ Plots the angle and magnitude of the 
           spatial averaged horizontal velocity (vx,vy) at every 
           level of the grid"""
        
        k = idxTime
        P = []
        T = []
        # Plot anlge over heigh
        if plotAngle:
            P += plot(self.zLevels, rad2deg(self.angleVxVy[k,:]),'-o')
            T += ["angle"]
            ylabel("angle in deg")
            if not angleRange is None:
                ylim(angleRange)
            
        xlabel("z in m")       
        if plotAngle and plotMag:
            twinx()

        # Plot mag over height
        if plotMag:
            P += plot(self.zLevels, self.magVxVy[k,:], 'o-r')
            T += ["mag"]
            ylabel("mag in m/s")
            if not magVeloRange is None:
                ylim(magVeloRange)

        # Set xrange and legend
        legend(P, T)
        


    def getAngleMagOverHeightImageSeq(self, prefixFileName,
                                      startTime, endTime,
                                      plotAngle = True,
                                      angleRange = None, plotMag = True,
                                      magVeloRange = None):
        """ Generate a list of images with the given prefixFileName 
          for a time intervall given by startTime and endTime in sec
          with the method plotAngleMagOverHeight. For the additional
          parameters look at this method."""

        if (magVeloRange is None) and plotMag:
            magVeloRangeL = self.getMagOverHeightRange(startTime, endTime)
        else:
            magVeloRangeL = magVeloRange
            
        for k in range(startTime, endTime):
            clf()
            title("time: %d s" % k)
            idx = searchsorted(self.time, k)
            self.plotAngleMagOverHeight(idx, plotAngle, angleRange,
                                        plotMag, magVeloRangeL) 
            savefig("%s_%06d.png" % (prefixFileName, k))


    def getMagOverHeightRange(self, startTime, endTime):
         """ Calculates the data range for the magnitude over height"""
         yMin = self.magVxVy[startTime:endTime, :].min()
         yMax = self.magVxVy[startTime:endTime, :].max()
         delta = yMax * 0.1
         yMin -= delta
         yMax += delta
         return yMin, yMax
         
            
    def plotVzEnergyOverHeight(self, idxTime, energyRange=None):
        """ Plots the vz-energy of the 
           spatial averaged horizontal velocity (vx,vy) at every 
           level of the grid"""
        
        k = idxTime

        plot(self.zLevels, self.sumVz[k,:], 'o-r', label="vz energy")
        ylabel("mag in m^2/s^2")
        if not energyRange is None:
                ylim(energyRange)

        xlabel("z in m")       
        legend()


    def getVzEnergyOverHeightRange(self, startTime, endTime):
        """ Get the maximum range of the plot for vz """
        yMin = self.sumVz[startTime:endTime, :].min()
        yMax = self.sumVz[startTime:endTime, :].max()
        delta = yMax * 0.1
        yMin -= delta
        yMax += delta
        return yMin, yMax
         
        
        
    def getVzEnergyOverHeightImageSeq(self, prefixFileName,
                                      startTime, endTime,
                                      energyRange=None):
        """ Generate a list of images with the given prefixFileName 
          for a time intervall given by startTime and endTime in sec
          with the method plotAngleMagOverHeight. For the additional
          parameters look at this method."""

        if energyRange is None:
            energyRangeL = self.getVzEnergyOverHeightRange(startTime, endTime)
        else:
            energyRangeL = energyRange
        
        for k in range(startTime, endTime):
            clf()
            title("time: %d s" % k)
            idx = searchsorted(self.time, k)
            self.plotVzEnergyOverHeight(idx, energyRange=energyRangeL) 
            savefig("%s_%06d.png" % (prefixFileName, k))

            
    def saveFileAnlgeMag(self, fileName):
        """ Saves the angle and magnitude of the 
           spatial averaged horizontal velocity (vx,vy) at every 
           level of the grid in the give file"""
        header = []
        columns = []

        # 1.) Column time
        header  += ["time in s"]
        columns += [self.time]
    
        # 2.) Go through all levels 
        N = len(self.zLevels)
        for i in range(N):
            angle = self.angleVxVy[:,i]
            mag   = self.magVxVy[:,i]
            s     = "at %0.3f m" % self.zLevels[i]

            header += ["angle " + s, "mag  " + s]
            columns += [array(angle), array(mag)]

        headerStr = reduce(lambda a,b: "%s, %s" % (a,b), header)
        A = column_stack(columns)
        savetxt(fileName, A, header = headerStr)


        
    def plotContourTimeZVz(self, nLevels=100):
        """ Generates a spatial-temporal plot of the 
            mag(vz) over height and time"""
        X,Y = np.meshgrid(self.time,self.zLevels)
        C = np.column_stack(self.sumVz)
        contourf(X,Y,C, levels=nLevels)
        xlabel("t in s")
        ylabel("z in m")
        colorbar()


        
    def plotContourTimeZMag(self, nLevels=100): 
        """ Generates a spatial-temporal plot of the 
            mag(vr) over height and time"""
        X,Y = np.meshgrid(self.time,self.zLevels)
        C = np.column_stack(self.magVxVy)
        contourf(X,Y,C, levels=nLevels)
        xlabel("t in s")
        ylabel("z in m")
        colorbar()

        
    def plotContourTimeAngle(self, nLevels=100): 
        """ Generates a spatial-temporal plot of the 
            angle(vr) over height and time"""
        X,Y = np.meshgrid(self.time,self.zLevels)
        C = np.column_stack(self.angleVxVy)
        contourf(X,Y,C, levels=nLevels)
        xlabel("t in s")
        ylabel("z in m")
        colorbar()

        
                
    def getChunksOfZLevels(self, numbersOfPartitions):
        """ Retuns a partion of the zLevels with the given number of 
        partitions."""
        N = self.zLevels.shape[0] - 2
        deltaN = N // numbersOfPartitions
        chunks = [(1+i*deltaN, 1+(i+1)*deltaN) \
                  for i in range(numbersOfPartitions)]
        chunks[-1] = (chunks[-1][0], max(chunks[-1][1], N))
        return chunks



    def calcAveragedVelocityForZLevelPartitions(self, numbersOfPartitions):
        """ Paritiones the meanVelociy and abs(vz) over the height over 
            the zLevels into partions. The number of paritions is given as 
            parameter """
        chunks = self.getChunksOfZLevels(numbersOfPartitions)
        
        # Partition vall 
        self.partitionVall = np.zeros((self.time.shape[0],
                                       numbersOfPartitions))
        # Partition vz
        self.partitionVz = np.zeros((self.time.shape[0], numbersOfPartitions))

        for i in range(numbersOfPartitions):
            self.partitionVall[:, i] = \
                self.meanVall[:,chunks[i][0]:chunks[i][1]].mean(axis=1) 
            self.partitionVz[:, i] = \
                self.absMeanVz[:,chunks[i][0]:chunks[i][1]].mean(axis=1) 
    

            
    def savePartitionsVz(self, fileName):
        """ Saves the partition data vz"""
        header = []

        for i in range(self.partitionVall.shape[1]):
            header  += ["Re abs(vz) partition %d" % i]
            
        headerStr = reduce(lambda a,b: "%s, %s" % (a,b), header)
        savetxt(fileName, self.partitionVz * self.ReFac, header = headerStr)

        
    def savePartitionsVall(self, fileName):
        """ Saves the partition data for vall"""
        header = []

        for i in range(self.partitionVall.shape[1]):
            header  += ["Re mag(v) partition %d" % i]
            
        headerStr = reduce(lambda a,b: "%s, %s" % (a,b), header)
        savetxt(fileName, self.partitionVall  * self.ReFac, header = headerStr)
