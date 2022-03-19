#import numpy as np
#import matplotlib.pyplot as plt
from numpy import *
from matplotlib.pyplot import *
import glob
import os
from multiprocessing import Pool 
from functools import reduce

class UDVAnalysisFelix:

    def __init__(self, fileName):
        self.data = np.loadtxt(fileName, skiprows=1)
        self.columns = ["t [s]",   
                        "u_lsc_top [m/s]",
                        "u_lsc_bot [m/s]",
                        "lsc_angle_top [deg]",
                        "lsc_angle_bot [deg]",
                        "t [t_turnover]",
                        "t normalized by <t_ff>",
                        "u_lsc_top [m/s]_normed on $<u_{ff}>_{t}$",
                        "u_lsc_bot [m/s]_normed on $<u_{ff}>_{t}$",
                        "u_lsc_top [m/s]_smoothed",
                        "u_lsc_bot [m/s]_smoothed",
                        "lsc_angle_top [deg]_smoothed",
                        "lsc_angle_bot [deg]_smoothed",
                        "u_lsc_top [m/s]_normed on $<u_{ff}>_{t}$_smoothed",
                        "u_lsc_bot [m/s]_normed on $<u_{ff}>_{t}$_smoothed"
                        ]

        
    def syncTime(self, timeOffset):
        self.data[:,0] -= timeOffset
        self.t0StartIndex = argwhere(self.data[:,0] >= 0)

        
    def plotULSC_Top(self):
        p = plot(self.data[:,0], self.data[:,1], label="u_lsc_topF")
        return p
        
    
    def plotULSC_Bottom(self):
        p = plot(self.data[:,0], self.data[:,2], label="u_lsc_botF")
        return p
    
    def plotAngleLSC_Top(self):
        p = plot(self.data[:,0], self.data[:,3], label="angle_lsc_topF")
        return p
    
    
    def plotAngleLSC_Bottom(self):
        p = plot(self.data[:,0], self.data[:,4], label="angle_lsc_botF")
        return p

    
    def plotTurnoverTime(self):
        p = plot(self.data[:,0], self.data[:,5], label="turnoverF")
        return p








    

class TemperatureAnalysisFelix:

    def __init__(self, fileName):
        self.data = np.loadtxt(fileName, skiprows=1)
        self.columns = ["t [s]",
                        "t normalized by <t_ff>",
                        "Ra",
                        "Pr",
                        "Cool_mean",
                        "Heat_mean",
                        "therm_boundary_layer_vert_1",
                        "therm_boundary_layer_vert_2",
                        "therm_boundary_layer_vert_3",
                        "t_diffusive [s]",
                        "t_diffusive [s]_smoothed",
                        "therm_boundary_layer_vert_1_smoothed",
                        "therm_boundary_layer_vert_2_smoothed",
                        "therm_boundary_layer_vert_3_smoothed"]

        
    def syncTime(self, timeOffset):
        self.data[:,0] -= timeOffset
        self.t0StartIndex = argwhere(self.data[:,0] >= 0)


    def plotRa(self):
        p = plot(self.data[:,0], self.data[:,2], label="Ra tempF")
        return p

    def plotPr(self):
        p = plot(self.data[:,0], self.data[:,3], label="Pr tempF")
        return p

    def plotThermBoundaryLayerVert1(self):
        p = plot(self.data[:,0], self.data[:,6],
                 label="thermalBoundary vert 1 F")
        return p

    def plotThermBoundaryLayerVert2(self):
        p = plot(self.data[:,0], self.data[:,7],
                 label="thermalBoundary vert 2 F")
        return p

    def plotThermBoundaryLayerVert3(self):
        p = plot(self.data[:,0], self.data[:,8],
                 label="thermalBoundary vert 3 F")
        return p

    def plotTDiffusive(self):
        plot(self.data[:,0], self.data[:,9], label="t-diffusive F")
        return p


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
            savefig("%s_%d.png" % (prefixFileName, k))
    
