class ExperimentConfigList:

    def __init__(self):
        """ Generates a database with the date as key"""
        self.database = {}

        
    def add(self, newExp):
        """ Adds one experiment """
        self.database[newExp.expDate] = newExp


    def listDates(self):
        """ Lists all saved dates """
        return list(self.database.keys())


    def getExperiment(self, theDate):
        """ Returns the experiment for the given date"""
        return self.database[theDate]


    def getListOfKeysForTemperatureDifference(self, deltaT):
        """ Returns a list of dates of experiments with the 
        same deltaT value"""
        theList = []
        for k in self.database:
            if self.database[k].deltaT == deltaT:
                theList += [k]

        return theList


    def getListOfKeysForSensorConfig(self, sensorConfigName):
        """ Returns a list of dates of experiments with the 
        same sensorConfig value"""
        theList = []
        for k in self.database:
            if self.database[k].sensorConfig == sensorConfigName:
                theList += [k]

        return theList


    def getListOfKeysForMagneticFieldConfig(self, magneticField):
        """ Returns a list of dates of experiments with the 
        same magnetic field configuration"""
        theList = []
        for k in self.database:
            if self.database[k].magneticField == magneticField:
                theList += [k]

        return theList
        

    def getListOfSensorConfigs(self):
        """ Returns all available sensor configs """
        theSet = set()
        for k in self.database:
            theSet.add(self.database[k].sensorConfig)

        return theSet

    
    def getListOfMagneticFieldConfigs(self):
        """ Returns all available sensor configs """
        theSet = set()
        for k in self.database:
            theSet.add(self.database[k].magneticField)

        return theSet

    
    def getListOfTemperatureDifference(self):
        """ Returns all available sensor temperature differences """
        theSet = set()
        for k in self.database:
            theSet.add(self.database[k].deltaT)

        return theSet


