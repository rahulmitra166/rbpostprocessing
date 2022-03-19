from experimentList import expConfList
from Experiment import Experiment
import sys

print("Usage: genVideo.py <YYYY_MM_DD> <Nr of partitions>")

if len(sys.argv) < 3:
    print("No date or nr of partitions for the experiment is given")
    sys.exit(-1)


# 1.) Select experiment
c      = expConfList.getExperiment(sys.argv[1])
theExp = Experiment(c)

# 2.) Load data
theExp.loadCIFTVelocity() 
theExp.loadTemp()
 
    
theExp.cift.calcAveragedVelocityForZLevelPartitions(int(sys.argv[2]))
theExp.cift.savePartitionsVz(theExp.getOutputPath("partitionVz.dat"))
theExp.cift.savePartitionsVall(theExp.getOutputPath("partitionVall.dat"))
