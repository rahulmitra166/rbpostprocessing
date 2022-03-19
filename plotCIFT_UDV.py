from matplotlib.pylab import *
from numpy import *
import os
import sys
import scipy.ndimage as ndimage
import glob

from DOPpy import *
from experimentList import expConfList

###############################################################################################################
print("Usage: plotCIFT_UDV.py <YYYY_MM_DD> <SensorToPlot>")

filtering = input("Filtering of UDV data? [g - Gauss/ t - Threshold/ n - No filtering] \n")

if(filtering == "t"):
  threshold = float(input("Set the threshold in mm/s \n"))
elif(filtering == "g"):
  sigma = float(input("Set the standard deviation of Gaussian kernel \n"))
  order = int(input("Set the order of Gaussian filter \n"))

depthLine = int(input("Enter the depth along which velocity will be plotted : "))

date = sys.argv[1]
channel_name = sys.argv[2]

channelID_list_39 = array([4,5,6,7,8,9,10])
channelID_list_59 = array([1,2,3,4,5,7,8,9,10])

channel_name_39 = array(["Top_0","34_0","Mid_0","14_0","Bot_0","Bot_45","Bot_135"])
channel_name_59 = array(["Mid_45","Mid_90","Mid_135","14_90","Bot_90","Top_90","34_90","Top_135","Top_45"])

UDOP = 0

if(channel_name in channel_name_39):
  print("Sensor is in UDOP39")
  UDOP = 39
elif(channel_name in channel_name_59):
  print("Sensor is in UDOP59")
  UDOP = 59
else:
  print("Error! Channel name not correct")
  sys.exit(-1)
###############################################################################################################

c = expConfList.getExperiment(date)

measurement_folder = "/mnt/archiv/Limcon/LargeRB_CIFT/01_measurements/"
output_folder = "/mnt/archiv/Limcon/LargeRB_CIFT/05_output/"
reko_folder = "/mnt/archiv/Limcon/LargeRB_CIFT/03_rekos/"

#--------------------------------------------------------------------------------------------------------------

UDV_data_folder =  measurement_folder+c.sensorConfig+"/"+date+"/"
plot_folder = output_folder+c.sensorConfig+"/"+date+"/cift_udvcopy/"
CIFT_data_folder = reko_folder+c.sensorConfig+"/"+date+"/CIFT_UDVcopy/"

###############################################################################################################
nTimeSteps = 30000
timeCIFT = linspace(0.0,nTimeSteps,nTimeSteps)

startO = c.cutInHours 
deltaT = c.deltaT

if (deltaT > 10):
  
  vmin = -24
  vmax = 24
  threshold = threshold * 1e-3
  lev = linspace(vmin,vmax,100)

if (deltaT < 10):

  vmin = -16
  vmax = 16
  threshold = threshold * 1e-3
  lev = linspace(vmin,vmax,100)


if(UDOP==39):
  
  print("Loading UDV files")
  fPattern39 = (UDV_data_folder + "%s_??_0.BDD" %(c.expDateShort))
  bdd39L = glob.glob(fPattern39)
  if len(bdd39L) == 0:
    print("No UDOP39 file found!")
  bdd39 = DOP(bdd39L[0])
  print("DOP39 channels:  ", bdd39.getChannels())
  
  channelID = take(channelID_list_39, where(channel_name_39==channel_name)[0])[0]
  
  time39 = bdd39.getTime(channelID)

  t1 = time39[1]
  t2 = time39[2]
  diff = t2-t1
  
  start39 = int((startO*3600)/diff)
  
  time = bdd39.getTime(channelID)
  data = bdd39.getChannelParam('velo', channelID)
  if(filtering == "t"):
    rows, cols = where(abs(data)>threshold)
    data[rows, cols] = nan
  elif(filtering == "g"):
    data = ndimage.gaussian_filter(data, sigma=sigma, order=order)
    data = data*1e3
  else:
    data = data
  
  time = time[start39:]
  data = data[start39:,:]
  data=data.T
  
  data = data*1e3
     
  time -= time[0]
    
  depth = bdd39.getDepth(channelID)
  Ndepth = len(depth)
      
  udv_velo_radial = load(CIFT_data_folder+"udv_velo_radial39_%d.npy" % channelID)
  udv_velo_radial = udv_velo_radial[:nTimeSteps,:]
  udv_velo_radial = udv_velo_radial.T
  udv_velo_radial = udv_velo_radial*1e3    
  
  figure(1)
  clf()
  gcf().set_size_inches([12,10])
      
  subplot(4,1,1)
  contourf(timeCIFT, depth, udv_velo_radial, cmap=cm.jet, levels=lev)
  
  if(filtering == "g"): 
    title(channel_name + " with gaussian filtering", fontsize=16)
  elif(filtering == "t"): 
    title(channel_name + " with threshold (%s mm/s) filtering" %(threshold*1e3), fontsize=16)
  else:
    title(channel_name + " without filtering", fontsize=16)
  ylabel("CIFT \n Depth (mm)")
  cb = colorbar()
  cb.set_ticks([lev.min(), (lev.min()/2), 0.0, (lev.max()/2), lev.max()])
      
  subplot(4,1,2)
  contourf(time, depth, data, cmap=cm.jet, levels=lev)
      
  ylabel("UDV \n Depth (mm)")
  xlim(0,nTimeSteps)
  cb = colorbar()
  cb.set_ticks([lev.min(), (lev.min()/2), 0.0, (lev.max()/2), lev.max()])
  
  subplot(4,1,3)
  plot(timeCIFT, udv_velo_radial[(Ndepth - depthLine),:], "r-")
  ylabel("CIFT \n Depth %dmm \n Velocity (mm/s)" %depthLine)
  xlim(0,nTimeSteps)
  
  
  subplot(4,1,4)
  plot(time, data[depthLine,:], "b-")
  ylabel("UDV \n Depth %dmm \n Velocity (mm/s)" %depthLine)
  xlabel("Time (s)")  
  xlim(0,nTimeSteps)
      
  savefig(plot_folder+"CIFT_UDV_%s.png" %channel_name)
      
  close(1)
  
if(UDOP==59):  
  
  print("Loading UDV files")
  fPattern59 = (UDV_data_folder + "%s_??_1.BDD" %(c.expDateShort))
  bdd59L = glob.glob(fPattern59)
  if len(bdd59L) == 0:
    print("No UDOP59 file found!")
  bdd59 = DOP(bdd59L[0])
  print("DOP59 channels:  ", bdd59.getChannels())
  
  channelID = take(channelID_list_59, where(channel_name_59==channel_name)[0])[0]
  time59 = bdd59.getTime(channelID)
  
  t1 = time59[1]
  t2 = time59[2]
  diff = t2-t1
  start59 = int((startO*3600)/diff)

  time = bdd59.getTime(channelID)
  data = bdd59.getChannelParam('velo', channelID)
  if(filtering == "t"):
    rows, cols = where(abs(data)>threshold)
    data[rows, cols] = nan
  elif(filtering == "g"):
    data = ndimage.gaussian_filter(data, sigma=sigma, order=order)
    data = data*1e3
  else:
    data = data

  time = time[start59:]
  data = data[start59:,:]
  data=data.T
  
  data = data*1e3   
  time -= time[0]
    
  depth = bdd59.getDepth(channelID)
  Ndepth = len(depth)
      
  udv_velo_radial = load(CIFT_data_folder+"udv_velo_radial59_%d.npy" % channelID)
  udv_velo_radial = udv_velo_radial[:nTimeSteps,:]
  udv_velo_radial = udv_velo_radial.T
  udv_velo_radial = udv_velo_radial*1e3    
  
  figure(1)
  clf()
  gcf().set_size_inches([12,10])
      
  subplot(4,1,1)
  contourf(timeCIFT, depth, udv_velo_radial, cmap=cm.jet, levels=lev)
      
  if(filtering == "g"): 
    title(channel_name + " with gaussian filtering", fontsize=16)
  elif(filtering == "t"): 
    title(channel_name + " with threshold (%s mm/s) filtering" %(threshold*1e3), fontsize=16)
  else:
    title(channel_name + " without filtering", fontsize=16)
  ylabel("CIFT \n Depth (mm)")
  cb = colorbar()
  cb.set_ticks([lev.min(), (lev.min()/2), 0.0, (lev.max()/2), lev.max()])
      
  subplot(4,1,2)
  contourf(time, depth, data, cmap=cm.jet, levels=lev)
      
  ylabel("UDV \n Depth (mm)")
  xlim(0,nTimeSteps)
  cb = colorbar()
  cb.set_ticks([lev.min(), (lev.min()/2), 0.0, (lev.max()/2), lev.max()])
  
  subplot(4,1,3)
  plot(timeCIFT, udv_velo_radial[(Ndepth - depthLine),:], "r-")
  ylabel("CIFT \n Depth %dmm \n Velocity (mm/s)" %depthLine)
  xlim(0,nTimeSteps)
  
  subplot(4,1,4)
  plot(time, data[depthLine,:], "b-")
  ylabel("UDV \n Depth %dmm \n Velocity (mm/s)" %depthLine)
  xlabel("Time (s)")  
  xlim(0,nTimeSteps)
      
  savefig(plot_folder+"CIFT_UDV_%s.png" %channel_name)
      
  close(1)
  