import numpy as np
from ExperimentConfigList import *
from ExperimentConfig import *


expConfList = ExperimentConfigList()


#########################################################################
#
# Experiments in August 2021
#
#########################################################################
S0 = "SensorConfig_00"
shiftElements = np.zeros(7,dtype=int)
shiftElements[0] = 1

expConfList.add(ExperimentConfig(expDate =  "2021_07_28", sensorConfig = S0,
                                 durationInH = 8, deltaT = 1.2, cutInHours = 2,
                                 magneticField = "Bz", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_07_29", sensorConfig = S0,
                                 durationInH = 8, deltaT = 1.2, cutInHours = 2,
                                 magneticField = "Bz", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_08_03", sensorConfig = S0,
                                 durationInH = 8, deltaT = 1.2, cutInHours = 3,
                                 magneticField = "Bz", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_08_04", sensorConfig = S0,
                                 durationInH = 8, deltaT = 12, cutInHours = 4,
                                 magneticField = "Bz", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_08_05", sensorConfig = S0,
                                 durationInH = 8, deltaT = 12, cutInHours = 4,
                                 magneticField = "Bz", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_08_09", sensorConfig = S0,
                                 durationInH = 8, deltaT = 12, cutInHours = 4,
                                 magneticField = "Bz", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_08_10", sensorConfig = S0,
                                 durationInH = 8, deltaT = 0.2, cutInHours = 4,
                                 magneticField = "Bz", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_08_11", sensorConfig = S0,
                                 durationInH = 8, deltaT = 0.2, cutInHours = 4,
                                 magneticField = "Bz", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_08_12", sensorConfig = S0,
                                 durationInH = 8, deltaT = 0.2, cutInHours = 4,
                                 magneticField = "Bz", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_08_16", sensorConfig = S0,
                                 durationInH = 8, deltaT = 1.2, cutInHours = 4,
                                 magneticField = "Bx", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_08_17", sensorConfig = S0,
                                 durationInH = 8, deltaT = 1.2, cutInHours = 4,
                                 magneticField = "Bx", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_08_18", sensorConfig = S0,
                                 durationInH = 8, deltaT = 1.2, cutInHours = 4,
                                 magneticField = "Bx", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_08_19", sensorConfig = S0,
                                 durationInH = 8, deltaT = 12, cutInHours = 4,
                                 magneticField = "Bx", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))
                                 
expConfList.add(ExperimentConfig(expDate =  "2021_08_31", sensorConfig = S0,
                                 durationInH = 8, deltaT = 12, cutInHours = 3,
                                 magneticField = "Bx", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_09_01", sensorConfig = S0,
                                 durationInH = 8, deltaT = 12, cutInHours = 4,
                                 magneticField = "Bz", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_09_02", sensorConfig = S0,
                                 durationInH = 8, deltaT = 0.2, cutInHours = 3,
                                 magneticField = "Bx", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_09_06", sensorConfig = S0,
                                 durationInH = 8, deltaT = 0.2, cutInHours = 4,
                                 magneticField = "Bx", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_09_07", sensorConfig = S0,
                                 durationInH = 8, deltaT = 0.2, cutInHours = 5,
                                 magneticField = "Bz", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_09_08", sensorConfig = S0,
                                 durationInH = 8, deltaT = 0.2, cutInHours = 4,
                                 magneticField = "Bx", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2021_09_09", sensorConfig = S0,
                                 durationInH = 8, deltaT = 12, cutInHours = 5,
                                 magneticField = "Bx", NusseltPT100 = False,
                                 flowRateCool = 160, flowRateHot = 160,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = shiftElements,
                                 comment = ""))
#########################################################################
#
# Experiments in January 2022
#
#########################################################################
S2 = "SensorConfig_02"
expConfList.add(ExperimentConfig(expDate =  "2022_01_20", sensorConfig = S2,
                                 durationInH = 8, deltaT = 1.2, cutInHours = 7,
                                 magneticField = "Bz", NusseltPT100 = True,
                                 flowRateCool = None, flowRateHot = None,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = None,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2022_01_24", sensorConfig = S2,
                                 durationInH = 8, deltaT = 12, cutInHours = 6,
                                 magneticField = "Bz", NusseltPT100 = True,
                                 flowRateCool = None, flowRateHot = None,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = None,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2022_02_01", sensorConfig = S2,
                                 durationInH = 8, deltaT = 2.6, cutInHours = 7,
                                 magneticField = "Bz", NusseltPT100 = True,
                                 flowRateCool = None, flowRateHot = None,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = None,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2022_02_02", sensorConfig = S2,
                                 durationInH = 8, deltaT = 5.35, cutInHours = 5,
                                 magneticField = "Bz", NusseltPT100 = True,
                                 flowRateCool = None, flowRateHot = None,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = None,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2022_02_07", sensorConfig = S2,
                                 durationInH = 8, deltaT = 7.0, cutInHours = 5,
                                 magneticField = "Bz", NusseltPT100 = True,
                                 flowRateCool = None, flowRateHot = None,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = None,
                                 comment = "flow starts already 21:00, lasts till 3:00"))

expConfList.add(ExperimentConfig(expDate =  "2022_02_08", sensorConfig = S2,
                                 durationInH = 8, deltaT = 9.0, cutInHours = 7,
                                 magneticField = "Bz", NusseltPT100 = True,
                                 flowRateCool = None, flowRateHot = None,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = None,
                                 comment = ""))
                                 
expConfList.add(ExperimentConfig(expDate =  "2022_02_09", sensorConfig = S2,
                                 durationInH = 8, deltaT = 21.5, cutInHours = 6,
                                 magneticField = "Bz", NusseltPT100 = True,
                                 flowRateCool = None, flowRateHot = None,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = None,
                                 comment = ""))

expConfList.add(ExperimentConfig(expDate =  "2022_02_14", sensorConfig = S2,
                                 durationInH = 8, deltaT = 14.18, cutInHours = 7.5,
                                 magneticField = "Bz", NusseltPT100 = True,
                                 flowRateCool = None, flowRateHot = None,
                                 nRings = 7, nSensorsAtRing = 6,
                                 shiftElements = None,
                                 comment = "experiment failed partially, reference sensor did not work"))

