from matplotlib.pyplot import *
from CIFTVelocityLargeRB import *
from ReadFelixAnalysis import *
import sys


pref = "1K" 
c = CIFTVelocityLargeRB("%s/Reko/"%pref, "%s/CIFTPoints/points_vol.dat"%pref) 

if pref == "1K":
    tempFelix = TemperatureAnalysisFelix("felixComparison/20210803_01_ra1e8_tc_timelines_cutdown.out")
    UDVFelix = UDVAnalysisFelix("felixComparison/20210803_01_udv_timelines_ra1e8_cutdown.out")
elif pref == "12K":
    tempFelix = TemperatureAnalysisFelix("felixComparison/20210809_01_ra1e9_tc_timelines_cutdown.out")
    UDVFelix = UDVAnalysisFelix("felixComparison/20210809_01_udv_timelines_ra1e9_cutdown_TW.out")
else:
    print("Not correct!")
    sys.exit(-1)
    
tempFelix.syncTime(3*3600)                                            
UDVFelix.syncTime(3*3600)                                             



######################################################################
#
# Comparison of angles
#
######################################################################
c.extractDirectionVectorAtZLevels(0.02)
angleBot = list(map(lambda s:rad2deg(s[1]),c.angleVxVyRaw))
angleTop = list(map(lambda s:rad2deg(s[23]),c.angleVxVyRaw))

clf()
plot(c.times, array(angleTop)+180, label="CIFT")
UDVFelix.plotAngleLSC_Top()
plot(UDVFelix.data[:,0], UDVFelix.data[:,11],label="angle smooth")
legend()
s = "top"
title(s)
savefig("%s_Uangle_%s_UDV.png" % (pref, s))
xlim(10000, 14000)
savefig("%s_Uangle_%s_UDV_zoom1.png" % (pref, s))
xlim(15000, 17000)
savefig("%s_Uangle_%s_UDV_zoom2.png" % (pref, s))


clf() 
plot(c.times, array(angleBot)+180, label="CIFT") 
UDVFelix.plotAngleLSC_Bottom() 
plot(UDVFelix.data[:,0], UDVFelix.data[:,12],label="angle smooth") 
legend() 
s = "bottom"
title(s)
savefig("%s_Uangle_%s_UDV.png" % (pref, s))
xlim(10000, 14000)
savefig("%s_Uangle_%s_UDV_zoom1.png" % (pref, s))
xlim(15000, 17000)
savefig("%s_Uangle_%s_UDV_zoom2.png" % (pref, s))


####################################################################
#
# Comparison of amplitudes
#
####################################################################
for j in [1,23]:
    clf()
    for k in [0.01, 0.02,0.05,0.1, 0.15, 0.2]:
        c.extractDirectionVectorAtZLevels(k)
        p = c.plotMagAtLevelOverTime(j)
        p[0].set_label("avg r: %0.3f" % k)
    if j == 1:
        UDVFelix.plotULSC_Bottom()
        s = "bottom"
    else:
        UDVFelix.plotULSC_Top()
        s = "top"
    legend()

    title(s)
    savefig("%s_Uavg_%s_UDV.png" % (pref, s))
    xlim(10000, 14000)
    savefig("%s_Uavg_%s_UDV_zoom1.png" % (pref, s))
    xlim(15000, 17000)
    savefig("%s_Uavg_%s_UDV_zoom2.png" % (pref, s))

    
    

########################################################################
#
# Quiver vxvy over height
#
#############äääääääääääääääääääääääääääääääääääääääääääääääääääääääääää
#clf()
#c.plotSpatiallyAveragedVelocity(True)
#twinx()
#a = tempFelix.plotThermBoundaryLayerVert1()[0]
#a.set_color("r") 

