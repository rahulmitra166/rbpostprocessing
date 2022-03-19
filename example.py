import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from experimentList import expConfList
from Experiment import Experiment

        
# 1.) Select experiment
c      = expConfList.getExperiment("2021_08_09")
theExp = Experiment(c)

# 2.) Load data
theExp.loadCIFTVelocity()
theExp.loadCIFTBField()
theExp.loadTemp()

showCIFTBField      = False
showCIFTVelocity    = False
showTemp            = True
showExperimentPlots = True


#############################################################################
#
# 3a.) Examples for bfield
#
#############################################################################
if showCIFTBField:
    # Figure 10
    theExp.bField.plotSpatialTemporalVertical([-200,200], 5, figNr = 10)
    title("Spatial-temporal plot of magnetic field at rings over height")
    
    # Figure 11
    theExp.bField.plotSpatialTemporalHorizontal([-300,300], figNr = 11)
    title("Spatial-temporal plot of magnetic field at rings over azimuth")

    figure(12); clf()
    title("Magnetic field over height at given time of 15000 s")
    theExp.bField.plotSensorsOverHeight(15000, [-300,300])

    # Figure 13
    theExp.bField.plotSinusFitAtRings(15000, [-150, 150],figNr = 13)
    title("Magnetic field with sin fit at rings")

    
    # Figure 14
    angle,amp = theExp.bField.getAmplAndAngleAtRings()
    fig, axes = plt.subplots(nrows=2, sharex=True, num=14, clear=True)
    fig.set_size_inches([16, max(8,4.4)])
    for i in range(amp.shape[1]):
        axes[0].plot(theExp.bField.times, angle[:,i], label="ring %d" % (i+1))
        axes[0].set_ylabel("angle in ra")
        axes[0].legend()
        axes[1].plot(theExp.bField.times, amp[:,i], label="ring %d" % (i+1))
        axes[1].set_ylabel("amplitude")
        axes[1].set_xlabel("t in s")
        axes[1].legend()


#############################################################################
#
# 3b.) Examples for CIFT velocity
#
#############################################################################
if showCIFTVelocity:
    # Figure 20 + 21
    figure(20); clf()
    title("Mean velocity over time in m/s")
    theExp.cift.plotSpatiallyAveragedVelocity(True)  # as Reynolds number
    figure(21); clf()
    title("Mean velocity over time as Reynolds number")
    theExp.cift.plotSpatiallyAveragedVelocity(False) # as velocity

    # Figure 22
    figure(22); clf()
    title("Max velocity over time as Reynolds number")
    theExp.cift.plotSpatiallyMaxVelocity(True)

    # Figure 23
    figure(23); clf()
    title("horizontal velocity at level 1 and t = 15000")
    theExp.cift.plotVectorVxVyAtLevel(15000, 1)

    # Figure 24
    figure(24); clf()
    title("horizontal velocity at top and bottom at t = 15000")
    theExp.cift.plotVectorVxVyTopBottom(15000)

    # Figure 25
    figure(25); clf()
    title("Angle of horizontal velocity at top and bottom")
    theExp.cift.plotAngleTopBottomOverTime()

    # Figure 26
    figure(26); clf()
    title("Angle of horizontal velocity at level 1 and 5")
    theExp.cift.plotAngleAtLevelOverTime(1)
    theExp.cift.plotAngleAtLevelOverTime(5)
    legend()

    # Figure 27
    figure(27); clf()
    title("Magnitude of horizontal velocity at top and bottom")
    theExp.cift.plotMagTopBottomOverTime()

    # Figure 28
    figure(28); clf()
    title("Magnitude of horizontal velocity at level 1 and 5")
    theExp.cift.plotMagAtLevelOverTime(1)
    theExp.cift.plotMagAtLevelOverTime(5)
    legend()

    # Figure 29
    figure(29); clf()
    title("Angle and magnitude of horizontal velocity over height at t = 15000")
    theExp.cift.plotAngleMagOverHeight(15000)

    # Figure 30
    figure(30); clf()
    title("Energy of vz over height at t = 15000")
    theExp.cift.plotVzEnergyOverHeight(15000)

    # Figure 31
    figure(31); clf()
    title("Spatial temporal plot of energy of vz over height and time")
    theExp.cift.plotContourTimeZVz()

    # Figure 32
    figure(32); clf()
    title("Spatial temporal plot of magnitude over height and time")
    theExp.cift.plotContourTimeZMag()

    # Figure 33
    figure(33); clf()
    title("Spatial temporal plot of angle over height and time")
    theExp.cift.plotContourTimeAngle()


    #theExp.cift.saveSpatiallyAveragedData("spatiallyAveragedData.dat")
    #theExp.cift.saveFileAnlgeMag("Anlge mag")



#############################################################################
#
# 3c.) Examples for temperature plots
#
#############################################################################
if showTemp:
    #theExp.temp.calcFreeFallTime()
    #theExp.temp.getRa()
    #theExp.temp.getPr()

    t =  theExp.temp.getTimeIntervall()
    print("Time intervall: %d s - %d s" % (t[0], t[1])),

    # Figure 40
    figure(40); clf() ; title("Interval with temperature difference")
    theExp.temp.checkIntervall()

    # Figure 41
    figure(41); clf() ; title("Nusselt Top")
    theExp.temp.plotNusselt(useCool = True)

    # Figure 42
    figure(42); clf() ; title("Nusselt Bottom")
    theExp.temp.plotNusselt(useCool = False)

    # Figure 43
    figure(43); clf() ; title("Temperature difference in bulk")
    theExp.temp.plotTemperatureDifferenceBulkOverTime()

    # Figure 44
    figure(44); clf() ; title("Temperature in bulk at 15000 s")
    theExp.temp.plotTemperatureDifferenceBulk(15000)
    
    # Figure 45
    figure(45); clf() ; title("Boundary layer 1")
    theExp.temp.plotThermBoundaryLayerVert1()

    # Figure 45
    figure(45); clf() ; title("Boundary layer 2")
    theExp.temp.plotThermBoundaryLayerVert2()

    # Figure 46
    figure(46); clf() ; title("Boundary layer 3")
    theExp.temp.plotThermBoundaryLayerVert3()

    # Figure 47
    figure(47); clf() ; title("Boundary layer 1")
    theExp.temp.plotThermBoundaryLayerNu()

    # Figure 48
    figure(48); clf() ; title("Plot diffusive time")
    theExp.temp.plotTDiffusive()

    # Save data
    #theExp.temp.writeTemperatureCorrectedData("example_tempCorrected.dat")
    #theExp.temp.saveSpatiallyAveragedData("example_spatiallyAveragedTemp.dat")

 

#############################################################################
#
# 3d.) Examples for temperature plots
#
#############################################################################
if showExperimentPlots:

    # Figure 46
    figure(46); clf() ; title("Re and temperature difference in the bulk")
    theExp.compareReAndTemperatureDiffInBulk()
