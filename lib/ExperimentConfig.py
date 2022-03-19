class ExperimentConfig:

    def __init__(self, expDate, sensorConfig, durationInH, deltaT, cutInHours,
                 magneticField, NusseltPT100, flowRateCool, flowRateHot,
                 nRings, nSensorsAtRing,
                 shiftElements = None, comment = ""):
        """ This class holds the information about one experiment 
        
        Arguments:
        expDate       -- Date of the experiment given as string "YYYY_MM_DD" 
        sensorConfig  -- Sensor config: directory name
        durationInH   -- duration of the experiment with temperature difference
                          in hours
        deltaT        -- Temperature difference
        cutInHours    -- Number of hours cut in the beginning of the experiment
        magneticField -- applied magnetic field
        NusseltPT100  -- use PT100 for Nusselt calculation at the top
        flowRateCool  -- Flowrate of the water of the cooling thermostat,  
                         None, if the flow rate was measured
        flowRateHot   -- Flowrate of the water of the heating thermostat,
                         None, if the flow rate was measured
        nRings        -- Number of sensor rings over height
        nSensorsAtRing -- 
        shiftElements  -- Shift azimutal postion, if the angle of the angle 
                          is not ascending for all sensors in the ring
                          (Needed for sine fit and plot over azimuth)
        comment        -- Any additional comment
        """
        
        # Experiment description
        self.expDate       = expDate

        self.sensorConfig  = sensorConfig
        self.durationInH   = durationInH
        self.deltaT        = deltaT
        self.cutInHours    = cutInHours
        self.comment       = comment
        
        # CIFT
        self.magneticField = magneticField
        self.radiusForVeloExtraction = 0.1
        self.shiftElements = shiftElements
        self.nRings = nRings # Sensor rings
        self.nSensorsAtRing = nSensorsAtRing

        
        # Temperature info
        self.NusseltPT100 = NusseltPT100
        self.flowRateCool = flowRateCool
        self.flowRateHot  = flowRateHot
        
        # Calculated parameters
        self.expDateShort  = self.expDate.replace("_","")


