from experimentList import expConfList
from Experiment import Experiment
import sys


print("Usage: comparisonReDeltaTNu.py <YYYY_MM_DD>")
if len(sys.argv) < 2:
    print("No date for the experiment is given")
    sys.exit(-1)
    
# 1.) Select experiment
c      = expConfList.getExperiment(sys.argv[1])
theExp = Experiment(c)

# 2.) Load data
theExp.loadCIFTVelocity() 
theExp.loadTemp()

# 3.) Plot special
theExp.compareReAndTemperatureDiffInBulk()

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    plt.show()
