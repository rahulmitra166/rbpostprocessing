from matplotlib.pyplot import *
from TemperatureProbesLargeRB import *
from CIFTVelocityLargeRB import *
from CIFTBFieldLargeRB import *
from ReadFelixAnalysis import *
import sys
import glob
from ExperimentConfig import *
from experimentList import expConfList
        

        
class Experiment:

    def __init__(self, expConfig):
        self.expConfig = expConfig
    
        # Directories
        self.measurementDir = "01_measurements"
        self.rekoDir        = "03_rekos"
        self.outputDir      = "05_output"


    def getMeasPath(self, postfix):
        """ Generates the path given by the experiment name for measurements"""
        return "%s/%s/%s/%s" % (self.measurementDir,self.expConfig.sensorConfig,
                                self.expConfig.expDate, postfix)

    def getRekoPath(self, postfix):
        """ Generates the path given by the experiment name for reko"""
        return "%s/%s/%s/%s" % (self.rekoDir,self.expConfig.sensorConfig,
                                self.expConfig.expDate, postfix)

    def getOutputPath(self, postfix):
        """ Generates the path given by the experiment name for output"""
        return "%s/%s/%s/%s" % (self.outputDir,self.expConfig.sensorConfig,
                                self.expConfig.expDate, postfix)
    

    def loadCIFTVelocity(self):
        """ Load CIFT reconstructed velocities """
        print("Load CIFT velocity")
        self.cift = CIFTVelocityLargeRB(self.getRekoPath("reko/"),
                                        self.getRekoPath("Reko_Setup")) 
        self.cift.extractDirectionVectorAtZLevels(    \
                self.expConfig.radiusForVeloExtraction)
   

    def loadCIFTBField(self):
        """ Load CIFT measured flow induced magnetic field """
        print("Load CIFT bfield")
        self.bField = CIFTBFieldLargeRB(self.getRekoPath("bfield/"),
                                        self.getRekoPath("Reko_Setup"),
                                        self.expConfig.nRings,
                                        self.expConfig.nSensorsAtRing,
                                        self.expConfig.shiftElements)


    def loadTemp(self):
        """ Load temperature files """
        print("Load temperature data")
        fPattern = self.getMeasPath(self.expConfig.expDateShort + "_??.dat")
        tempDataFileL = glob.glob(fPattern)

        if len(tempDataFileL) == 0:
            print("No temperature file found!")
            return

        tempDataFile = tempDataFileL[0]
        self.temp = TemperatureProbesLargeRB(tempDataFile)
        self.temp.calcNusselt(self.expConfig.NusseltPT100,
                              self.expConfig.flowRateCool,
                              self.expConfig.flowRateHot)
        
        self.temp.syncTime(self.expConfig.cutInHours * 3600)
        self.timeInterval = self.temp.getTimeIntervall()
                              

    def loadUDV(self):
        print("Load UDV data")
        


    def saveData(self):
        # Save data for further analysis
        self.cift.saveSpatiallyAveragedData( \
                        self.getOutputPath("velocitySpatialAverage.dat"))
        self.cift.saveFileAnlgeMag(self.getOutputPath("velocityAngleMag.dat"))
        self.temp.saveSpatiallyAveragedData( \
                self.getOutputPath("tempSpatialAveraged.dat"))
        self.temp.writeTemperatureCorrectedData( \
                self.getOutputPath("tempCorrected.dat"))


    def genImagesForVideo(self, timeIntervall, startTime = None,
                          endTime = None):
        """ Generates a list of images for the video """
        startTime = startTime if not startTime is None else self.timeInterval[0]
        endTime = endTime if not endTime is None else self.timeInterval[1]

        c = self.cift
        
        print("Generate image sequence Re plot")
        c.getSpatiallyMeanMaxVelocityImageSeq(self.getOutputPath("movie/meanV"),
                                              startTime, endTime, timeIntervall,
                                              plotMaxVelocity=False)
        
        print("Generate image sequence for angle over height")
        c.getAngleMagOverHeightImageSeq(self.getOutputPath("movie/angle"),
                                startTime, endTime, plotAngle=True,
                                plotMag=False)
        
        print("Generate image sequence for amplitude over height")
        c.getAngleMagOverHeightImageSeq(self.getOutputPath("movie/mag"),
                                        startTime, endTime, plotAngle=False,
                                        plotMag=True)
        
        print("Generate image sequence for vz-energy ")
        c.getVzEnergyOverHeightImageSeq(self.getOutputPath("movie/energyvz"),
                                        startTime, endTime)
        
        print("Generate image sequence for horizontal velocity top/bottom ")
        c.getVectorVxVyTopBottomImageSeq(self.getOutputPath("movie/vecVxVz"),
                                         startTime, endTime)

                         
        print("Generate image sequence for thermal boundary thickness")
        self.temp.getVzThermalBoundaryImageSeq( \
                            self.getOutputPath("movie/thermalB"),
                            startTime, endTime,
                            timeIntervall)

        print("Generate image sequence for bulk temperature")
        self.temp.plotTemperatureDifferenceBulkImageSeq( \
                            self.getOutputPath("movie/bulkT"),
                            startTime, endTime,
                            yRange = None)
        
        print("Generate image sequence for Nusselt number")
        self.temp.getNusseltImageSeq(self.getOutputPath("movie/Nusselt"),
                                     startTime, endTime,
                                     timeIntervall, False)

        print("Write temperature VTK files")
        self.temp.writeVTKFileSeq(self.getOutputPath("vtk/temp"), 0, 30000)




    #######################################################################
    #
    # Combine all images to one file
    #
    #######################################################################
    def combineImages(self, three3DViewPrefix,
                      startTime = None, endTime = None):
        """ Combines the images for each time step into one image for
            generating the movie """
        
        startTime = startTime if not startTime is None else self.timeInterval[0]
        endTime = endTime if not endTime is None else self.timeInterval[1]
        for i in range(startTime, endTime):
            b11 = self.getOutputPath("movie/meanV_%06d.png" % i)
            b21 = self.getOutputPath("movie/Nusselt_%06d.png"%  i)
            #b31 = self.getOutputPath("movie/thermalB_%06d.png" % i) 
            b31 = self.getOutputPath("movie/bulkT_%06d.png" % i) 
            b12 = self.getOutputPath("movie/mag_%06d.png" % i)
            b22 = self.getOutputPath("movie/angle_%06d.png" % i)
            b32 = self.getOutputPath("movie/energyvz_%06d.png" % i)
            
            zout = self.getOutputPath("movie/zout.png")

            os.system(("montage %s %s %s %s %s %s "%(b11,b12,b21,b22,b31,b32))+\
                      "-tile 2x3 -mode Concatenate %s" % zout )
            
            c1 = self.getOutputPath("movie/vecVxVz_%06d.png" %  i)

            if not three3DViewPrefix is None:
                c2 = self.getOutputPath("movie/%s%04d.png" %
                                        (three3DViewPrefix,i))
            else:
                c2 = c1 # Do nothing else
        
            newFile =  self.getOutputPath("movie/movie_%06d.png" %  i)
            os.system("montage %s %s %s -mode Concatenate -trim "%
                      (zout, c1,c2)+\
                      "-scale 1280x1440+0+0 -geometry +1+1 %s" % newFile)


            
    #######################################################################
    #
    # Combine all images to one file
    #
    #######################################################################
    def genMovieFromImages(self, startTime = None, endTime = None,
                        frameRate = 25):
        startTime = startTime if not startTime is None else self.timeInterval[0]
        endTime = endTime if not endTime is None else self.timeInterval[1]

        outputMovie = self.getOutputPath("movie/movie.mp4")
        imageList   = self.getOutputPath("movie/movie_*.png")

        cmd = ("ffmpeg -framerate %d -pattern_type glob " % frameRate) + \
            ("-i \"%s\" " % imageList) + \
            "-c:v libx264 -profile:v high " + \
            " -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" " + \
            " -y -crf 20 -pix_fmt " + \
            "yuv420p %s" % outputMovie
        print(cmd)
        os.system(cmd)


    def compareReAndTemperatureDiffInBulk(self):
        """ Generates a plot for comparison of Re with thermal
            difference in the bulk"""
        clf()
        gcf().set_size_inches([17.97,  4.8 ])
        self.cift.plotSpatiallyAveragedVelocity(True)
        L = gca().get_lines()
        twinx()
        self.temp.plotTemperatureDifferenceBulkOverTime()
        L += gca().get_lines()
        L[0].set_color("r")
        L[0].set_label("Re")
        legend(L, list(map(lambda s:s.get_label(), L)))

 
