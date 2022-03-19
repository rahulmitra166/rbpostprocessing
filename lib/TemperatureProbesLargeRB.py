import numpy as np
import matplotlib.pyplot as plt
from numpy import *
from matplotlib.pyplot import *
import glob
import os
from multiprocessing import Pool 
from functools import reduce
from scipy import interpolate
from configparser import ConfigParser






class TemperatureProbesLargeRB:
    """ This class reads the recorded temperature data and 
        does the analyis
    """
    
    def __init__(self, fileNameTemperature):
        self.fileNameTemperature = fileNameTemperature

        # 1.) Load temperature file
        D = np.loadtxt(fileNameTemperature,skiprows=100, delimiter=";")
        self.data = D
        # Get inf-Columns
        self.infCols = list(argwhere(np.isinf(self.data).sum(axis=0) > 0)[:,0])

        self._loadConfig()

        # 1a.) Init sensor configuration
        self._initSensorConfig()
        self._applyCalibration()
        self._removeInfValuesOfAllSensorRings()
    
        # 2.) get time
        self.t = D[:,1]

        # 2a.) get temperature different top bottom
        self.deltaTGaInSnSimp = D[:,self.sensorsBottomCuPlate].mean(axis=1) - \
                                D[:,self.sensorsTopCuPlate].mean(axis=1)


        # 3.) Init some parameters
        self.nusseltCool = np.zeros(self.t.shape[0])
        self.nusseltHot  = np.zeros(self.t.shape[0])
        self._calcMeanTemperatures()
        self._calcLocalThermalBoundary()


        

    def _loadConfig(self):
        fName =  self.fileNameTemperature.replace(".dat", ".ini")
        if not os.path.isfile(fName):
            print("File %s not found!" % fName)
            return
            
        conf = ConfigParser()
        conf.read(fName)

        self.channelNoToSensName = {}
        for i in conf.sections(): 
            if i.startswith("Temperature"):
                channel = int(conf[i]['Channel'])
                name    = conf[i]['Name'].replace("\"","")
                used    = bool(conf[i]['Used'])
                if used:
                    self.channelNoToSensName[channel] = name  
 


    def _initSensorConfig(self):
        self._initSensorNames()
        self._setSensorPositionsAndGroups()

        
    def _initSensorNames(self):
        # ColumnID -> Sensor Name
        self.sensorNames = \
            ["Vol Flux Cool","time",
             "Vol Flux Heat", "Cooling bath",  "Heating bath", "Heating Power",
             "PT100_1_res", "PT100_2_res", "PT100_1_T", "PT100_2_T",
             "Bot01_rect", "Mid01_rect", "Top01_rect",
             "Top02","Top03","Top04","Top05","Top06","Top07","Top08",
             "Mid10","Mid11","Mid12","Mid13","Mid14","Mid15","Mid16",
             "Bot02","Bot03","Bot04","Bot05","Bot06","Bot07","Bot08",
             "Vert03_layer","Vert03_bulk","Vert01_layer","Vert01_bulk",
             "Vert02_layer","Vert02_bulk",
             "Vert07","Vert04","Vert06","Vert05",
             "Hall01",
             "Cool01","Cool02","Cool03","Cool04","Cool05","Cool06","Cool07",
             "Cool08",
             "Top09","Top10","Top11","Top12","Top13","Top14","Top15","Top16",
             "Mid02","Mid03","Mid04","Mid05","Mid06","Mid07","Mid08","Mid09",
             "Bot09","Bot10","Bot11","Bot12","Bot13","Bot14","Bot15","Bot16",
             "Cool_in","Cool_out","Heat_in","Heat_out",
             "Hall02",
             "Heat01","Heat02","Heat03","Heat04","Heat05","Heat06","Heat07",
             "Heat08"]


        # Sensor Name -> ColumnID
        self.sensorMapToIndex = {}
        for i,cnt in zip(self.sensorNames, range(len(self.sensorNames))):
            self.sensorMapToIndex[i] = cnt

                            
    def _setSensorPositionsAndGroups(self):
        """ Set sensor information """
        self.sensorsTopCuPlate   = [45,46,47,48,49,50,51,52]  # 8 sonden
        self.sensorsBottomCuPlate = [82,83,84,85,86,86,88,89]  # 8 sonden
        self.sensorsVertical     = [41,43,42,40]  # Vert04:Vert07 (layer)

        self.sensorRingTop       = [13,14,15,16,17,18,19,
                                    53,54,55,56,57,58,59,60] # start with top2
        self.sensorRingMid       = [61,62,63,64,65,66,67,68,
                                    20,21,22,23,24,25,26] # start with mid2
        self.sensorRingBot       = [27,28,29,30,31,32,33,
                                    69,70,71,72,73,74,75,76] # start with bot2
        self.hall                = [44,81]

        # Sensor position in coordinate system of Felix:
        # Sensors top plate
        zCu = 0.3255
        self.posSensorsTopCuPlate = \
            np.array([[ 0.092261019848815,  0.07673268023767, zCu],
                      [ 0.010980194239608,  0.119496591308959, zCu],
                      [-0.07673268023767,   0.092261019848815, zCu],
                      [-0.119496591308959,  0.010980194239608, zCu],
                      [-0.092261019848815, -0.07673268023767, zCu],
                      [-0.010980194239608, -0.119496591308959, zCu],
                      [ 0.07673268023767,  -0.092261019848815, zCu],
                      [ 0.119496591308959, -0.010980194239608, zCu]])

        # Sensors bottom plate
        self.posSensorsBottomCuPlate = self.posSensorsTopCuPlate.copy()
        self.posSensorsBottomCuPlate[:,2] = 	-0.325
        
        # Sensors vertical
        self.posSensorsVertical = \
             np.array([[0, 0, 0.32],
                       [0, 0, 0.32],
                       [0, 0, 0.32],
                       [0, 0, 0.32]])

        
        # Sensor ring top positions
        zt = 0.16
        self.posSensorRingTop = \
            np.array([[ 0.133035137968407,  0.088891237283136, zt],
                      [ 0.088891237283136,  0.133035137968407, zt],
                      [ 0.031214451522581,  0.156925644864517, zt],
                      [-0.031214451522581,  0.156925644864517, zt],
                      [-0.088891237283136,  0.133035137968407, zt],
                      [-0.133035137968407,  0.088891237283136, zt],
                      [-0.156925644864517,  0.031214451522581, zt],
                      [-0.156925644864517, -0.031214451522581, zt],
                      [-0.133035137968407, -0.088891237283136, zt],
                      [-0.088891237283137, -0.133035137968407, zt],
                      [-0.031214451522581, -0.156925644864517, zt],
                      [ 0.031214451522581, -0.156925644864517, zt],
                      [ 0.088891237283136, -0.133035137968407, zt],
                      [ 0.133035137968407, -0.088891237283136, zt],
                      [ 0.156925644864517, -0.031214451522581, zt]])
        
        # Sensor ring mid positions
        self.posSensorRingMid = self.posSensorRingTop.copy()
        self.posSensorRingMid[:,2] = 0

        # Sensor ring bot positions
        self.posSensorRingBot = self.posSensorRingTop.copy()
        self.posSensorRingBot[:,2] = -0.16

        # Transfer to CIFT Coordinate system
        def toCIFT(D):
            R = np.zeros_like(D)
            angle = deg2rad(180)
            M = np.array([[cos(angle), -sin(angle), 0],
                          [sin(angle),  cos(angle), 0],
                          [0,           0,          1]])
            for i in range(D.shape[0]):
                R[i,:] = np.dot(M, D[i,:])
            return R
        
        self.posSensorsTopCuPlate    = toCIFT(self.posSensorsTopCuPlate)
        self.posSensorsBottomCuPlate = toCIFT(self.posSensorsBottomCuPlate)
        self.posSensorsVertical      = toCIFT(self.posSensorsVertical)
        self.posSensorRingTop        = toCIFT(self.posSensorRingTop)
        self.posSensorRingMid        = toCIFT(self.posSensorRingMid)
        self.posSensorRingBot        = toCIFT(self.posSensorRingBot)


        
    def _removeInfValues(self, refMap):
        """ Removes the inf values by interplolation for the 
            given subset (refMap). refMap is a sensor ring. 
            The circular nature of the data is taken into account
        """
        cols    = set(self.infCols).intersection(set(refMap))
        infList = []
        okList  = []
        for i,k in zip(refMap, range(len(refMap))):
            if i in cols:
                infList += [k]
            else:
                okList += [k]
                
        if len(infList) > 0:
            # Cirular interpolation!
            x = list(map(lambda x:x-len(refMap), okList)) + \
                okList + \
                list(map(lambda x:x+len(refMap), okList)) 
                
            for k in range(self.data.shape[0]):
                y = list(self.data[k, list(map(lambda j: refMap[j], okList))])
                y = y + y + y
                f = interpolate.interp1d(x, y)
                for i in range(len(infList)):
                        self.data[k, refMap[infList[i]]] = f(infList[i])

                        
        
    def _removeInfValuesOfAllSensorRings(self):
        """ Removes inf information from the temperature data"""
        self._removeInfValues(self.sensorRingTop)
        self._removeInfValues(self.sensorRingMid)
        self._removeInfValues(self.sensorRingBot)
                
 
            
    def _applyCalibration(self):
        """ Apply the calibration to all sensors"""
        self._polynomCorrection()
        self._offsetCorrection()


        
    def _polynomCorrection(self):
        """ Does the polynomial correction """
        fName = "04001_01_20191101_tc_polynomials.pol"
        rootDir = os.environ.get("LIMCON_POST_ROOT")
        if rootDir is None:
            print("Environment variable not set: LIMCON_POST_ROOT")
            print("%s not loaded" % fName)
            return

        # Load polynomial correction and save it in the map
        # chanell -> polynom
        fullName = rootDir + "/tempCalibration/" + fName
        polyStr =  loadtxt(fullName, skiprows=1,dtype='str')
        polyStr = polyStr.transpose()[1:,:]
        self.polyMap = {}
        for i in range(polyStr.shape[0]):
            self.polyMap[self.channelNoToSensName[int(polyStr[i,0])]] = \
                list(map(float, polyStr[i,1:]))

        # Helper function to calculate the 3rd order polynomial
        polyF = lambda x,k: x**3 * k[0] + x**2 * k[1] + x * k[2] + k[3]
        
        # correction of the temprature data    
        for k in self.polyMap:
            dataCol = self.data[:, self.sensorMapToIndex[k]]
            mask    = np.isfinite(dataCol) # ignore inf values
            
            self.data[mask, self.sensorMapToIndex[k]] = \
                polyF(dataCol[mask], self.polyMap[k])
            


        
    def _offsetCorrection(self):
        """ Does the offset correction according to Felix """
        fName = "20191209_02_offset.csv"
        rootDir = os.environ.get("LIMCON_POST_ROOT")
        if rootDir is None:
            print("Environment variable not set: LIMCON_POST_ROOT")
            print("%s not loaded" % fName)
            return

        fullName = rootDir + "/tempCalibration/" + fName
        offsetStr = loadtxt(fullName, skiprows=1,dtype='str')
        self.offsetMap = {}
        for i in range(offsetStr.shape[0]):
            self.offsetMap[offsetStr[i,0]] = float(offsetStr[i,1])
            
        for k in self.offsetMap:
            self.data[:, self.sensorMapToIndex[k]] += self.offsetMap[k]
        

            
    def _calcMeanTemperatures(self):
        """ Calcuates the mean temperature at the sensor rings """
        self.meanRingTop = self.data[:, self.sensorRingTop].mean(axis=1)
        self.meanRingMid = self.data[:, self.sensorRingMid].mean(axis=1)
        self.meanRingBot = self.data[:, self.sensorRingBot].mean(axis=1)
        self.meanTopCuPlate = self.data[:, self.sensorsTopCuPlate].mean(axis=1)
        self.meanBottomCuPlate=self.data[:,
                                         self.sensorsBottomCuPlate].mean(axis=1)

        
    def _calcLocalThermalBoundary(self):
        """ Caclulate the thermal boundaries at sensors 
        Vert01_layer, Vert02_layer, Vert03_layer"""

        t_0mm = self.data[:,self.sensorsVertical].mean(axis=1)
        t_mean =(self.data[:,self.sensorsBottomCuPlate].mean(axis=1) + \
                 self.data[:,self.sensorsTopCuPlate].mean(axis=1)) / 2.0
        
        t_3mm = self.data[:,self.sensorMapToIndex["Vert01_layer"]]
        self.boundaryLayer1 = 0.003*(t_0mm-t_mean)/(t_0mm-t_3mm) * 1e3
        
        t_3mm = self.data[:,self.sensorMapToIndex["Vert02_layer"]]
        self.boundaryLayer2 = 0.003*(t_0mm-t_mean)/(t_0mm-t_3mm)  * 1e3
        
        t_3mm = self.data[:,self.sensorMapToIndex["Vert02_layer"]]
        self.boundaryLayer3 = 0.003*(t_0mm-t_mean)/(t_0mm-t_3mm) * 1e3

        self.boundaryLayerNu = np.zeros_like(self.boundaryLayer1)  * 1e3
        
        
    def calcFreeFallTime(self):
        pass

    def getRa(self):
        pass

    def getPr(self):
        pass
        
        
    def calcNusselt(self, highPrecisionTempProbes,
                    constFlowRateCool = None, constFlowRateHot = None):
        """ Calculates the time dependent Nusselt number of the 
            temperature measurements. 
            
            Arguments:
        highPrecisionTempProbes -- if true, use high precision PT 100, otherwise
                                   the standard thermocouples
        constFlowRateCool -- give the flow rate for cooling, if not measured
        constFlowRateHot -- give the flow rate for heating, if not measured
        """
        D = self.data
        
        # 0a) Phsical constants for GaInSn at 20 degree C
        g     = 9.81    # m/s^2
        self.g = g
        
        # GaInSn
        alphaGaInSn = 1.24e-4 # 1/K (thermal expansion coefficient, or beta)
        self.alphaGaInSn = alphaGaInSn
        rhoGaInSn   = 6350    #  kg/m^3
        cpGaInSn    = 366     # J/(kg K)
        nuGaInSn    = 3.44e-7 # m^2/s
        lambdaGaInSn= 24.0    # W/(m K)

        # Water at 20 degree C
        alphaH2O = 69e-6    # 1/K (thermal expansion coefficient, or beta)
        rhoH2O   = 997      #  kg/m^3
        cpH2O    = 4184     # J/(kg K)
        nuH2O    = 1.0e-6   # m^2/s
        lambdaH2O= 0.143    # W/(m K)

        # 0b) vessel:
        d    = 0.320
        H    = 2*d
        area =  1.0/4.0 * pi * d**2
        self.d = d
        self.H = H
        
        
        
        flowrC = D[:,0] if constFlowRateCool is None else constFlowRateCool
        self.coolFlow =  flowrC  / 1795.0 * 1e-3 # l/s - m^3/s
        
        if highPrecisionTempProbes:
            self.coolIn   = D[:,8] 
            self.coolOut  = D[:,9] 
        else:
            self.coolIn   = D[:,77]
            self.coolOut  = D[:,78]

        flowrH = D[:,2] if constFlowRateHot is None else constFlowRateHot
        self.hotFlow = flowrH  / 915.0 * 1e-3 # l/s - m/s
        
        if highPrecisionTempProbes:
            self.hotIn   = D[:,79]  # No special columns!
            self.hotOut  = D[:,80]  # 
        else:
            self.hotIn   = D[:,79]
            self.hotOut  = D[:,80]

        # 3a.) Time dependent temperature difference in GaInSn
        self.deltaTGaInSnComp = D[:,self.sensorsBottomCuPlate].mean(axis=1) + \
                                D[:,self.sensorsTopCuPlate].mean(axis=1)-2.0 * \
                                D[:,self.sensorsVertical].mean(axis=1)


        self.deltaTGaInSn = self.deltaTGaInSnComp

        # 3b.) Time dependent temperature difference in water
        self.deltaTH2OCool  = self.coolIn - self.coolOut 
        self.deltaTH2OHot  = self.hotIn - self.hotOut 

        # 4.) Heat conduction
        self.heatConductionGaInSn = lambdaGaInSn * self.deltaTGaInSn * (area/H)
        self.heatH2OCool          = cpH2O * rhoH2O * self.deltaTH2OCool * \
                                    self.coolFlow
        self.heatH2OHot          = cpH2O * rhoH2O * self.deltaTH2OHot * \
                                    self.hotFlow

        self.nusseltCool = self.heatH2OCool / self.heatConductionGaInSn
        self.nusseltHot = self.heatH2OHot / self.heatConductionGaInSn

        # 5.) Extract parameters
        self.calcDiffusiveTime()
        self.extractParameters()

        # 6.) Calc boundary layer based on Nusselt number
        self.boundaryLayerNu[self.startTidx:self.endTidx] = \
            self.H / (2.0 * self.nusseltCool[self.startTidx:self.endTidx]) \
            * 1e3


        
    def extractParameters(self, threashold = 0.2):
        """ Extract the time intervall with a temperature difference"""
        self.startTidx = 0
        self.endTidx   = self.t.shape[0] - 1
        
        threshold = self.deltaTGaInSnSimp.max() * 0.8
        idx = np.argwhere(np.diff(self.deltaTGaInSnSimp > threshold,
                                  prepend=False))[:,0]
        if idx.shape[0] < 2:
            print("\n WARNING: Intervall not found!\n")
            return
        else:
            self.startTidx = idx[0]
            self.endTidx   = idx[-1]

        # Calculate mean temperature
        self.meanDeltaT = self.deltaTGaInSnSimp[self.startTidx:
                                                self.endTidx].mean()
        self.t_ff = np.sqrt(self.H / self.g * self.alphaGaInSn *
                            self.meanDeltaT)

        
            
    def checkIntervall(self):
        """ Plots the temperature difference and marks the 
            automaticly detected start and end points"""
        plt.clf()
        plt.plot(self.t, self.deltaTGaInSnComp)
        plt.axvline(self.t[self.startTidx], c="red")
        plt.axvline(self.t[self.endTidx], c="red")
        plt.grid()


        
    def syncTime(self, timeOffset):
        """ Synchronise the time by subtracting the given time offset
            from the time """
        self.t -= timeOffset
        iList = argwhere(self.t >= 0)
        if len(iList) > 0:
            self.t0StartIndex = iList[0]
        else:
            self.t0StartIndex = None
            

    def getTimeIntervall(self):
        """ Returns the interesting interval based, where the 
            temperature difference is switched on"""
        return [int(round(self.t[self.startTidx])),
                int(round(self.t[self.endTidx]))]
   
        
    def plotNusselt(self, useCool = True):
        """ Plots the nusselt number over time"""
        if useCool:
            nusselt = self.nusseltCool
        else:
            nusselt = self.nusseltHot
        if not self.startTidx is None:
            plt.plot(self.t[self.startTidx:self.endTidx],
                     nusselt[self.startTidx:self.endTidx],
                     label="Nusselt")
        else:
            plt.plot(self.t, nusselt, label="Nusselt")
        
        
    def getNusseltImageSeq(self, prefixFileName,
                                      startTime, endTime,
                                      deltaT, useCool = False):
        """ Generate a list of images with the given prefixFileName 
          for a time intervall given by startTime and endTime in sec
          with the method plotNusselt. For the additional
          parameters look at this method."""

        
        clf()
        self.plotNusselt(useCool)
        legend()
        verticalLine = axvline(startTime)
        title("nusselt " + "top" if useCool else "bottom")
        for k in range(startTime, endTime):
            title("time: %d s" % k)
            verticalLine.set_data([k, k], [0,1])
            xlim(k - deltaT, k + deltaT)
            draw()
            savefig("%s_%06d.png" % (prefixFileName, k))



    def plotTemperatureDifferenceBulkOverTime(self):
        """ Plots the temperature difference in the bulk over
            time """
        DeltaBot = self.meanRingBot - self.meanRingMid 
        DeltaTop = self.meanRingMid - self.meanRingTop
        DeltaAll = self.meanRingBot - self.meanRingTop
        
        plot(self.t, DeltaBot, label="Delta bottom")
        plot(self.t, DeltaTop, label="Delta top")
        plot(self.t, DeltaAll, label="Delta full")



    def plotTemperatureDifferenceBulk(self, t):
        """ Plots the spatially averaged temperature at the 
            three rings of the vessel for a given time t"""
        k = searchsorted(self.t, t)  
        plot(['top', 'mid', 'bot'], 
             [self.meanRingTop[k], self.meanRingMid[k],
              self.meanRingBot[k]], '-o')
        
        

    def plotTemperatureDifferenceBulkImageSeq(self, prefixFileName,
                                              startTime, endTime,
                                              yRange = None):
        """ Generate a list of images with the given prefixFileName 
          for a time intervall given by startTime and endTime in sec
          with the method plotTemperatureDifferenceBulk. For the additional
          parameters look at this method."""

        if yRange is None:
            startIDX = searchsorted(self.t, startTime)
            endIDX   = searchsorted(self.t, endTime)
            yRangeL = [self.meanRingTop[startIDX:endIDX].min(),
                       self.meanRingBot[startIDX:endIDX].max()]
        else:
            yRangeL = xRange
            
        for k in range(startTime, endTime):
            clf()
            title("time: %d s" % k)
            self.plotTemperatureDifferenceBulk(k)
            ylim(yRangeL)
            savefig("%s_%06d.png" % (prefixFileName, k))


    def calcDiffusiveTime(self):
        # k = waermleitfaehigkeit.gainsn_geberth(heat_mean,cool_mean)/(
        #                   form.cp.gainsn_geberth(heat_mean,cool_mean)*
        #    form.rho.gainsn_geberth(heat_mean,cool_mean))
        #= lambda/(cp*rho)
        #t_diff = boundary_layer_nu**2/k
        self.diffusiveTime = np.zeros_like(self.t)
        
        # rho = Dichte GaInSn, kannste auch konstante nehmen, ich habs halt
        # abhängig von der Mittleren Temperatur im Bulk ausrechnen lassen (kann
        # ich alles auch in die Bibliothek einfügen).
        # cp = Spezifische Wärmekapazität GaInSn
        # waermleitfaehigkeit = Wärmeleitfähigkeit lambda von GaInSn
        

    def plotThermBoundaryLayerVert1(self):
        """ Plot boundary layer 1 """
        p = plot(self.t[self.startTidx:self.endTidx],
                 self.boundaryLayer1[self.startTidx:self.endTidx],
                 label="thermalBoundary vert 1")
        return p

    def plotThermBoundaryLayerVert2(self):
        """ Plot boundary layer 2 """
        p = plot(self.t[self.startTidx:self.endTidx],
                 self.boundaryLayer2[self.startTidx:self.endTidx],
                 label="thermalBoundary vert 2")
        return p

    def plotThermBoundaryLayerVert3(self):
        """ Plot boundary layer 3 """
        p = plot(self.t[self.startTidx:self.endTidx],
                 self.boundaryLayer3[self.startTidx:self.endTidx],
                 label="thermalBoundary vert 3")
        return p

    
    def plotThermBoundaryLayerNu(self):
        """ Plot boundary layer based on Nusselt number """
        p = plot(self.t[self.startTidx:self.endTidx],
                 self.boundaryLayerNu[self.startTidx:self.endTidx],
                 label="thermalBoundary Nusselt")
        return p

    
    def plotTDiffusive(self):
        """ Plots the diffusive time """
        pass
        #p = plot(self.data[:,0], self.data[:,9], label="t-diffusive F")
        #return p


    
    def getVzThermalBoundaryImageSeq(self, prefixFileName,
                                      startTime, endTime,
                                      deltaT):
        """ Generate a list of images with the given prefixFileName 
          for a time intervall given by startTime and endTime in sec
          with the method plotThermBoundaryLayerVert1. For the additional
          parameters look at this method."""
        clf()
        self.plotThermBoundaryLayerVert1()
        self.plotThermBoundaryLayerVert2()
        self.plotThermBoundaryLayerVert3()
        legend()
        verticalLine = axvline(startTime)

        for k in range(startTime, endTime):
            title("time: %d s" % k)
            verticalLine.set_data([k,k], [0,1])
            xlim(k - deltaT, k + deltaT)
            draw()
            savefig("%s_%06d.png" % (prefixFileName, k))
    
        
    def writeTemperatureCorrectedData(self, fileName):
       """ Writes out the temperature compensated data """
       savetxt(fileName, self.data)

       
    def saveSpatiallyAveragedData(self, fileName):
       """ Saves the spatial averaged data """
       A = np.column_stack((self.t,
                            self.meanTopCuPlate, self.meanBottomCuPlate,
                            self.meanRingTop, self.meanRingMid,
                            self.meanRingBot,
                            self.nusseltCool, self.nusseltHot,
                            self.boundaryLayer1,
                            self.boundaryLayer2,
                            self.boundaryLayer3,
                            self.boundaryLayerNu,
                            self.diffusiveTime))
       
       headerList = ["time in s",
                     "mean temp at top cu plate in degC",
                     "mean temp at bottom cu plate in degC",
                     "mean temp top ring in degC",
                     "mean temp mid ring in degC",
                     "mean temp bot ring in degC",
                     "Nussel top","Nusselt bottom",
                     "Local boundary layer 1 in mm",
                     "Local boundary layer 2 in mm",
                     "Local boundary layer 3 in mm",
                     "Diffusive time in s"]
       
       headerStr = reduce(lambda a,b:a + ", " + b, headerList)
       savetxt(fileName, A, header=headerStr)


                   

    def writeVTKFile(self, fileNamePrefix, t, topBottom=False):
        if topBottom:
            N = len(self.posSensorsTopCuPlate) + \
                len(self.posSensorsBottomCuPlate) + \
                len(self.posSensorsVertical)
        else:
            N = len(self.posSensorRingTop) + \
                len(self.posSensorRingMid) + \
                len(self.posSensorRingBot)
            
        k = searchsorted(self.t, t)
        
        f = open("%s_%06d.vtk" % (fileNamePrefix, t), "w")
        f.write("# vtk DataFile Version 2.0\n")
        f.write("Temperature: %d\n" % t)
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write("POINTS %d float\n" % N)
        if topBottom:
            for i in range(len(self.posSensorsTopCuPlate)):
                f.write("%e %e %e\n" % tuple(self.posSensorsTopCuPlate[i,:]))
            for i in range(len(self.posSensorsBottomCuPlate)):
                f.write("%e %e %e\n" % tuple(self.posSensorsBottomCuPlate[i,:]))
            for i in range(len(self.posSensorsVertical)):
                f.write("%e %e %e\n" % tuple(self.posSensorsVertical[i,:]))
        else:
            for i in range(len(self.posSensorRingTop)):
                f.write("%e %e %e\n" % tuple(self.posSensorRingTop[i,:]))
            for i in range(len(self.posSensorRingMid)):
                f.write("%e %e %e\n" % tuple(self.posSensorRingMid[i,:]))
            for i in range(len(self.posSensorRingBot)):
                f.write("%e %e %e\n" % tuple(self.posSensorRingBot[i,:]))

        f.write("POINT_DATA %d\n" % N)
        f.write("SCALARS temperature float 1\n")
        f.write("LOOKUP_TABLE default\n")
        
        #remInf = lambda s:0 if np.isinf(s) else s
        remInf = lambda s:s
        if topBottom:
            for i in range(len(self.posSensorsTopCuPlate)):
                f.write("%e\n"%remInf(self.data[k, self.sensorsTopCuPlate[i]]))
            for i in range(len(self.posSensorsBottomCuPlate)):
                f.write("%e\n"%remInf(self.data[k,
                                                self.sensorsBottomCuPlate[i]]))
            for i in range(len(self.posSensorsVertical)):
                f.write("%e\n" % remInf(self.data[k, self.sensorsVertical[i]]))
        else:
            for i in range(len(self.posSensorRingTop)):
                f.write("%e\n" % remInf(self.data[k, self.sensorRingTop[i]]))
            for i in range(len(self.posSensorRingMid)):
                f.write("%e\n" % remInf(self.data[k, self.sensorRingMid[i]]))
            for i in range(len(self.posSensorRingBot)):
                f.write("%e\n" % remInf(self.data[k, self.sensorRingBot[i]]))
        f.close()


    def writeVTKFileSeq(self, fileNamePrefix, tStart, tEnd, topBottom=False):
        for t in range(tStart, tEnd):
            self.writeVTKFile(fileNamePrefix, t, topBottom)
