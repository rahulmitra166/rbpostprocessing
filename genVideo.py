from experimentList import expConfList
from Experiment import Experiment
import sys


print("Usage: genVideo.py <YYYY_MM_DD>")
if len(sys.argv) < 2:
    print("No date for the experiment is given")
    sys.exit(-1)
        
# 1.) Select experiment
c      = expConfList.getExperiment(sys.argv[1])
theExp = Experiment(c)

# 2.) Load data
theExp.loadCIFTVelocity() 
theExp.loadTemp()

# Generate movie
startT, endT = theExp.temp.getTimeIntervall()
#startT = 15000
#endT   = 15100


theExp.genImagesForVideo(200, startT, endT)
theExp.combineImages(None, startT, endT)
theExp.genMovieFromImages(startT, endT)

