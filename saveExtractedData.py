from experimentList import expConfList
from Experiment import Experiment
import sys


print("Usage: extractData.py <YYYY_MM_DD>")
if len(sys.argv) < 2:
    print("No date for the experiment is given")
    sys.exit(-1)
        
# 1.) Select experiment
c      = expConfList.getExperiment(sys.argv[1])
theExp = Experiment(c)

# 2.) Load data
theExp.loadCIFTVelocity() 
theExp.loadTemp()

# Save data
print("Save data")
theExp.saveData() 

print("Write temperature VTK files")
theExp.temp.writeVTKFileSeq(theExp.getOutputPath("vtk/temp"), 0, 30000)
