# state file generated using paraview version 5.10.0

# uncomment the following three lines to ensure this script works in future versions
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

##########################################################################
#
# Thomas
#
##########################################################################
dateStr          = "2021_08_03"
veloField        = 'VeloReko'
scaleFactorVelo  = 8

def getVTKFileList(localPrefix):
   import glob
   import re
   import os

   # 1.) Search files with absolute path
   globStr = os.getcwd() + "/" + localPrefix
   fileList = glob.glob(globStr)

   # 2.) Sort files with ending numbers
   fileList.sort(key=lambda var:[int(x) if x.isdigit() else x 
                                 for x in re.findall(r'[^0-9]|[0-9]+', var)]) 

   return fileList 

fileNamesTemp   = getVTKFileList("05_output/%s/vtk/temp*.vtk" % dateStr) 
fileNamesBField = getVTKFileList("03_rekos/%s/Bfield/*.vtk" % dateStr) 
fileNamesVelo   = getVTKFileList("03_rekos/%s/Reko/*.vtk" % dateStr) 




##########################################################################

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [635, 782]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-1.3511135626988366, -0.3712094920291416, 0.6718928381627969]
renderView1.CameraFocalPoint = [-0.022515359460240823, 0.05464905191290771, -0.009742887872573658]
renderView1.CameraViewUp = [0.42780914695730865, 0.10269481730111009, 0.898016207136125]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.40189155537291377
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(635, 782)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
temp___ = LegacyVTKReader(registrationName='temp___*', FileNames=fileNamesTemp)

# create a new 'Legacy VTK Reader'
bfield_sensor_velo_ = LegacyVTKReader(registrationName='bfield_sensor_velo_*', FileNames=fileNamesBField)

# create a new 'Glyph'
sensorBox = Glyph(registrationName='sensorBox', Input=bfield_sensor_velo_,
    GlyphType='Box')
sensorBox.OrientationArray = ['POINTS', 'sensorNormals']
sensorBox.ScaleArray = ['POINTS', 'No scale array']
sensorBox.ScaleFactor = 0.02
sensorBox.GlyphTransform = 'Transform2'
sensorBox.GlyphMode = 'All Points'

# create a new 'Glyph'
bfieldVektor = Glyph(registrationName='bfieldVektor', Input=bfield_sensor_velo_,
    GlyphType='Arrow')
bfieldVektor.OrientationArray = ['POINTS', 'sensorNormals']
bfieldVektor.ScaleArray = ['POINTS', 'bfield_sensor']
bfieldVektor.ScaleFactor = 0.002
bfieldVektor.GlyphTransform = 'Transform2'
bfieldVektor.GlyphMode = 'All Points'

# create a new 'Glyph'
tempBoxes = Glyph(registrationName='TempBoxes', Input=temp___,
    GlyphType='Box')
tempBoxes.OrientationArray = ['POINTS', 'No orientation array']
tempBoxes.ScaleArray = ['POINTS', 'No scale array']
tempBoxes.ScaleFactor = 0.03
tempBoxes.GlyphTransform = 'Transform2'
tempBoxes.GlyphMode = 'All Points'

# create a new 'Legacy VTK Reader'
veloReko_ = LegacyVTKReader(registrationName='veloReko_*', FileNames=fileNamesVelo)

# create a new 'Vortex Cores'
vortexCores1 = VortexCores(registrationName='VortexCores1', Input=veloReko_)
vortexCores1.VectorField = ['POINTS', veloField]

# create a new 'Tube'
tube5 = Tube(registrationName='Tube5', Input=vortexCores1)
tube5.Scalars = ['POINTS', 'delta-criterion']
tube5.Vectors = ['POINTS', 'acceleration']
tube5.Radius = 0.005866659879684448

# create a new 'Glyph'
veloVectors = Glyph(registrationName='veloVectors', Input=veloReko_,
    GlyphType='Arrow')
veloVectors.OrientationArray = ['POINTS', veloField]
veloVectors.ScaleArray = ['POINTS', veloField]
veloVectors.ScaleFactor = scaleFactorVelo
veloVectors.GlyphTransform = 'Transform2'
veloVectors.MaximumNumberOfSamplePoints = 4000
veloVectors.Seed = 800

# create a new 'Stream Tracer'
streamTracer2 = StreamTracer(registrationName='StreamTracer2', Input=veloReko_,
    SeedType='Line')
streamTracer2.Vectors = ['POINTS', veloField]
streamTracer2.MaximumStreamlineLength = 0.6399999856948853

# init the 'Line' selected for 'SeedType'
streamTracer2.SeedType.Point1 = [0.14, 0.0, -0.3199999928474426]
streamTracer2.SeedType.Point2 = [0.14, 0.0, 0.3199999928474426]
streamTracer2.SeedType.Resolution = 20

# create a new 'Stream Tracer'
streamTracer3 = StreamTracer(registrationName='StreamTracer3', Input=veloReko_,
    SeedType='Line')
streamTracer3.Vectors = ['POINTS', veloField]
streamTracer3.MaximumStreamlineLength = 0.6399999856948853

# init the 'Line' selected for 'SeedType'
streamTracer3.SeedType.Point1 = [0.0, -0.14, -0.3199999928474426]
streamTracer3.SeedType.Point2 = [0.0, -0.14, 0.3199999928474426]
streamTracer3.SeedType.Resolution = 20

# create a new 'Tube'
tube3 = Tube(registrationName='Tube3', Input=streamTracer3)
tube3.Scalars = ['POINTS', 'AngularVelocity']
tube3.Vectors = ['POINTS', 'Normals']
tube3.Radius = 0.0031999999284744265

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=veloReko_,
    SeedType='Line')
streamTracer1.Vectors = ['POINTS', veloField]
streamTracer1.MaximumStreamlineLength = 0.6399999856948853

# init the 'Line' selected for 'SeedType'
streamTracer1.SeedType.Point1 = [-0.14, 0.0, -0.3199999928474426]
streamTracer1.SeedType.Point2 = [-0.14, 0.0, 0.3199999928474426]
streamTracer1.SeedType.Resolution = 20

# create a new 'Tube'
tube1 = Tube(registrationName='Tube1', Input=streamTracer1)
tube1.Scalars = ['POINTS', 'AngularVelocity']
tube1.Vectors = ['POINTS', 'Normals']
tube1.Radius = 0.0031999999284744265

# create a new 'Stream Tracer'
streamTracer4 = StreamTracer(registrationName='StreamTracer4', Input=veloReko_,
    SeedType='Line')
streamTracer4.Vectors = ['POINTS', veloField]
streamTracer4.MaximumStreamlineLength = 0.6399999856948853

# init the 'Line' selected for 'SeedType'
streamTracer4.SeedType.Point1 = [0.0, 0.14, -0.3199999928474426]
streamTracer4.SeedType.Point2 = [0.0, 0.14, 0.3199999928474426]
streamTracer4.SeedType.Resolution = 20

# create a new 'Tube'
tube4 = Tube(registrationName='Tube4', Input=streamTracer4)
tube4.Scalars = ['POINTS', 'AngularVelocity']
tube4.Vectors = ['POINTS', 'Normals']
tube4.Radius = 0.0031999999284744265

# create a new 'Tube'
tube2 = Tube(registrationName='Tube2', Input=streamTracer2)
tube2.Scalars = ['POINTS', 'AngularVelocity']
tube2.Vectors = ['POINTS', 'Normals']
tube2.Radius = 0.0031999999284744265

# create a new 'Annotate Time'
annotateTime1 = AnnotateTime(registrationName='AnnotateTime1')
annotateTime1.Format = 'Time: {time:0.0f} s'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from veloReko_
veloReko_Display = Show(veloReko_, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
veloReko_Display.Representation = 'Surface'
veloReko_Display.ColorArrayName = [None, '']
veloReko_Display.Opacity = 0.26
veloReko_Display.SelectTCoordArray = 'None'
veloReko_Display.SelectNormalArray = 'None'
veloReko_Display.SelectTangentArray = 'None'
veloReko_Display.OSPRayScaleArray = veloField
veloReko_Display.OSPRayScaleFunction = 'PiecewiseFunction'
veloReko_Display.SelectOrientationVectors = veloField
veloReko_Display.ScaleFactor = 0.06399999856948853
veloReko_Display.SelectScaleArray = 'None'
veloReko_Display.GlyphType = 'Arrow'
veloReko_Display.GlyphTableIndexArray = 'None'
veloReko_Display.GaussianRadius = 0.0031999999284744265
veloReko_Display.SetScaleArray = ['POINTS', veloField]
veloReko_Display.ScaleTransferFunction = 'PiecewiseFunction'
veloReko_Display.OpacityArray = ['POINTS', veloField]
veloReko_Display.OpacityTransferFunction = 'PiecewiseFunction'
veloReko_Display.DataAxesGrid = 'GridAxesRepresentation'
veloReko_Display.PolarAxes = 'PolarAxesRepresentation'
veloReko_Display.ScalarOpacityUnitDistance = 0.04372717621481763
veloReko_Display.OpacityArrayName = ['POINTS', veloField]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
veloReko_Display.ScaleTransferFunction.Points = [-0.000679604010656476, 0.0, 0.5, 0.0, 0.0007535369950346649, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
veloReko_Display.OpacityTransferFunction.Points = [-0.000679604010656476, 0.0, 0.5, 0.0, 0.0007535369950346649, 1.0, 0.5, 0.0]

# show data from bfield_sensor_velo_
bfield_sensor_velo_Display = Show(bfield_sensor_velo_, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'bfield_sensor'
bfield_sensorLUT = GetColorTransferFunction('bfield_sensor')
bfield_sensorLUT.RGBPoints = [-57.95941162109375, 0.02, 0.3813, 0.9981, -54.60271780831473, 0.02000006, 0.424267768, 0.96906969, -51.246023995535715, 0.02, 0.467233763, 0.940033043, -47.889330182756694, 0.02, 0.5102, 0.911, -44.53263636997768, 0.02000006, 0.546401494, 0.872669438, -41.17594255719867, 0.02, 0.582600362, 0.83433295, -37.81924874441964, 0.02, 0.6188, 0.796, -34.462554931640625, 0.02000006, 0.652535156, 0.749802434, -31.105861118861608, 0.02, 0.686267004, 0.703599538, -27.74916730608259, 0.02, 0.72, 0.6574, -24.392473493303577, 0.02000006, 0.757035456, 0.603735359, -21.03577968052455, 0.02, 0.794067037, 0.55006613, -17.679085867745535, 0.02, 0.8311, 0.4964, -14.322392054966514, 0.021354336738172372, 0.8645368555261631, 0.4285579460761159, -10.9656982421875, 0.023312914349117714, 0.897999359924484, 0.36073871343115577, -7.609004429408479, 0.015976108242848862, 0.9310479513349017, 0.2925631815088092, -4.252310616629465, 0.27421074700988196, 0.952562960995083, 0.15356836602739213, -0.8956168038504444, 0.4933546281681699, 0.9619038625309482, 0.11119493614749336, 2.4610770089285694, 0.6439, 0.9773, 0.0469, 5.81777082170759, 0.762401813, 0.984669591, 0.034600153, 9.174464634486597, 0.880901185, 0.992033407, 0.022299877, 12.531158447265625, 0.9995285432627147, 0.9995193706781492, 0.0134884641450013, 15.887852260044653, 0.999402998, 0.955036376, 0.079066628, 19.244546072823667, 0.9994, 0.910666223, 0.148134024, 22.60123988560268, 0.9994, 0.8663, 0.2172, 25.957933698381694, 0.999269665, 0.818035981, 0.217200652, 29.314627511160722, 0.999133332, 0.769766184, 0.2172, 32.671321323939736, 0.999, 0.7215, 0.2172, 36.02801513671875, 0.99913633, 0.673435546, 0.217200652, 39.384708949497764, 0.999266668, 0.625366186, 0.2172, 42.74140276227679, 0.9994, 0.5773, 0.2172, 46.098096575055806, 0.999402998, 0.521068455, 0.217200652, 49.45479038783482, 0.9994, 0.464832771, 0.2172, 52.81148420061383, 0.9994, 0.4086, 0.2172, 56.16817801339286, 0.9947599917687346, 0.33177297300202935, 0.2112309638520206, 59.524871826171875, 0.9867129505479589, 0.2595183410914934, 0.19012239549291934, 62.88156563895089, 0.9912458875646419, 0.14799417507952672, 0.21078892136920357, 66.2382594517299, 0.949903037, 0.116867171, 0.252900603, 69.59495326450893, 0.903199533, 0.078432949, 0.291800389, 72.95164707728796, 0.8565, 0.04, 0.3307, 76.30834089006694, 0.798902627, 0.04333345, 0.358434298, 79.66503470284599, 0.741299424, 0.0466667, 0.386166944, 83.021728515625, 0.6837, 0.05, 0.4139]
bfield_sensorLUT.ColorSpace = 'RGB'
bfield_sensorLUT.NanColor = [1.0, 0.0, 0.0]
bfield_sensorLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
bfield_sensor_velo_Display.Representation = 'Surface'
bfield_sensor_velo_Display.ColorArrayName = ['POINTS', 'bfield_sensor']
bfield_sensor_velo_Display.LookupTable = bfield_sensorLUT
bfield_sensor_velo_Display.SelectTCoordArray = 'None'
bfield_sensor_velo_Display.SelectNormalArray = 'None'
bfield_sensor_velo_Display.SelectTangentArray = 'None'
bfield_sensor_velo_Display.OSPRayScaleArray = 'bfield_sensor'
bfield_sensor_velo_Display.OSPRayScaleFunction = 'PiecewiseFunction'
bfield_sensor_velo_Display.SelectOrientationVectors = 'sensorNormals'
bfield_sensor_velo_Display.ScaleFactor = 0.0479999989271164
bfield_sensor_velo_Display.SelectScaleArray = 'bfield_sensor'
bfield_sensor_velo_Display.GlyphType = 'Arrow'
bfield_sensor_velo_Display.GlyphTableIndexArray = 'bfield_sensor'
bfield_sensor_velo_Display.GaussianRadius = 0.0023999999463558195
bfield_sensor_velo_Display.SetScaleArray = ['POINTS', 'bfield_sensor']
bfield_sensor_velo_Display.ScaleTransferFunction = 'PiecewiseFunction'
bfield_sensor_velo_Display.OpacityArray = ['POINTS', 'bfield_sensor']
bfield_sensor_velo_Display.OpacityTransferFunction = 'PiecewiseFunction'
bfield_sensor_velo_Display.DataAxesGrid = 'GridAxesRepresentation'
bfield_sensor_velo_Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
bfield_sensor_velo_Display.ScaleTransferFunction.Points = [-2.7361268997192383, 0.0, 0.5, 0.0, 1.4507410526275635, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
bfield_sensor_velo_Display.OpacityTransferFunction.Points = [-2.7361268997192383, 0.0, 0.5, 0.0, 1.4507410526275635, 1.0, 0.5, 0.0]

# show data from temp___
temp___Display = Show(temp___, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'temperature'
temperatureLUT = GetColorTransferFunction('temperature')
temperatureLUT.RGBPoints = [22.27692985534668, 0.231373, 0.298039, 0.752941, 22.390409469604492, 0.865003, 0.865003, 0.865003, 22.503889083862305, 0.705882, 0.0156863, 0.14902]
temperatureLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
temp___Display.Representation = 'Surface'
temp___Display.ColorArrayName = ['POINTS', 'temperature']
temp___Display.LookupTable = temperatureLUT
temp___Display.SelectTCoordArray = 'None'
temp___Display.SelectNormalArray = 'None'
temp___Display.SelectTangentArray = 'None'
temp___Display.OSPRayScaleArray = 'temperature'
temp___Display.OSPRayScaleFunction = 'PiecewiseFunction'
temp___Display.SelectOrientationVectors = 'None'
temp___Display.ScaleFactor = 0.031999999284744264
temp___Display.SelectScaleArray = 'temperature'
temp___Display.GlyphType = 'Arrow'
temp___Display.GlyphTableIndexArray = 'temperature'
temp___Display.GaussianRadius = 0.0015999999642372132
temp___Display.SetScaleArray = ['POINTS', 'temperature']
temp___Display.ScaleTransferFunction = 'PiecewiseFunction'
temp___Display.OpacityArray = ['POINTS', 'temperature']
temp___Display.OpacityTransferFunction = 'PiecewiseFunction'
temp___Display.DataAxesGrid = 'GridAxesRepresentation'
temp___Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
temp___Display.ScaleTransferFunction.Points = [22.458049774169922, 0.0, 0.5, 0.0, 22.53519058227539, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
temp___Display.OpacityTransferFunction.Points = [22.458049774169922, 0.0, 0.5, 0.0, 22.53519058227539, 1.0, 0.5, 0.0]

# show data from tempBoxes
tempBoxesDisplay = Show(tempBoxes, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tempBoxesDisplay.Representation = 'Surface'
tempBoxesDisplay.ColorArrayName = ['POINTS', 'temperature']
tempBoxesDisplay.LookupTable = temperatureLUT
tempBoxesDisplay.SelectTCoordArray = 'None'
tempBoxesDisplay.SelectNormalArray = 'Normals'
tempBoxesDisplay.SelectTangentArray = 'None'
tempBoxesDisplay.OSPRayScaleArray = 'temperature'
tempBoxesDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
tempBoxesDisplay.SelectOrientationVectors = 'None'
tempBoxesDisplay.ScaleFactor = 0.03519999980926514
tempBoxesDisplay.SelectScaleArray = 'temperature'
tempBoxesDisplay.GlyphType = 'Arrow'
tempBoxesDisplay.GlyphTableIndexArray = 'temperature'
tempBoxesDisplay.GaussianRadius = 0.001759999990463257
tempBoxesDisplay.SetScaleArray = ['POINTS', 'temperature']
tempBoxesDisplay.ScaleTransferFunction = 'PiecewiseFunction'
tempBoxesDisplay.OpacityArray = ['POINTS', 'temperature']
tempBoxesDisplay.OpacityTransferFunction = 'PiecewiseFunction'
tempBoxesDisplay.DataAxesGrid = 'GridAxesRepresentation'
tempBoxesDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tempBoxesDisplay.ScaleTransferFunction.Points = [22.458049774169922, 0.0, 0.5, 0.0, 22.529619216918945, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tempBoxesDisplay.OpacityTransferFunction.Points = [22.458049774169922, 0.0, 0.5, 0.0, 22.529619216918945, 1.0, 0.5, 0.0]

# show data from veloVectors
veloVectorsDisplay = Show(veloVectors, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for veloField
veloRekoLUT = GetColorTransferFunction(veloField)
veloRekoLUT.RGBPoints = [-0.01940820924937725, 0.02, 0.3813, 0.9981, -0.018627270923129147, 0.02000006, 0.424267768, 0.96906969, -0.017846332596881048, 0.02, 0.467233763, 0.940033043, -0.01706539427063295, 0.02, 0.5102, 0.911, -0.016284455944384848, 0.02000006, 0.546401494, 0.872669438, -0.015503517618136745, 0.02, 0.582600362, 0.83433295, -0.014722579291888645, 0.02, 0.6188, 0.796, -0.013941640965640545, 0.02000006, 0.652535156, 0.749802434, -0.013160702639392445, 0.02, 0.686267004, 0.703599538, -0.012379764313144345, 0.02, 0.72, 0.6574, -0.011598825986896244, 0.02000006, 0.757035456, 0.603735359, -0.010817887660648139, 0.02, 0.794067037, 0.55006613, -0.01003694933440004, 0.02, 0.8311, 0.4964, -0.009256011008151939, 0.021354336738172372, 0.8645368555261631, 0.4285579460761159, -0.008475072681903839, 0.023312914349117714, 0.897999359924484, 0.36073871343115577, -0.0076941343556557375, 0.015976108242848862, 0.9310479513349017, 0.2925631815088092, -0.006913196029407638, 0.27421074700988196, 0.952562960995083, 0.15356836602739213, -0.006132257703159536, 0.4933546281681699, 0.9619038625309482, 0.11119493614749336, -0.005351319376911436, 0.6439, 0.9773, 0.0469, -0.004570381050663335, 0.762401813, 0.984669591, 0.034600153, -0.0037894427244152366, 0.880901185, 0.992033407, 0.022299877, -0.0030085043981671333, 0.9995285432627147, 0.9995193706781492, 0.0134884641450013, -0.00222756607191903, 0.999402998, 0.955036376, 0.079066628, -0.0014466277456709302, 0.9994, 0.910666223, 0.148134024, -0.0006656894194228304, 0.9994, 0.8663, 0.2172, 0.00011524890682526948, 0.999269665, 0.818035981, 0.217200652, 0.0008961872330733728, 0.999133332, 0.769766184, 0.2172, 0.0016771255593214726, 0.999, 0.7215, 0.2172, 0.0024580638855695724, 0.99913633, 0.673435546, 0.217200652, 0.0032390022118176723, 0.999266668, 0.625366186, 0.2172, 0.004019940538065776, 0.9994, 0.5773, 0.2172, 0.004800878864313875, 0.999402998, 0.521068455, 0.217200652, 0.005581817190561975, 0.9994, 0.464832771, 0.2172, 0.006362755516810075, 0.9994, 0.4086, 0.2172, 0.007143693843058178, 0.9947599917687346, 0.33177297300202935, 0.2112309638520206, 0.007924632169306278, 0.9867129505479589, 0.2595183410914934, 0.19012239549291934, 0.008705570495554378, 0.9912458875646419, 0.14799417507952672, 0.21078892136920357, 0.009486508821802478, 0.949903037, 0.116867171, 0.252900603, 0.010267447148050581, 0.903199533, 0.078432949, 0.291800389, 0.011048385474298681, 0.8565, 0.04, 0.3307, 0.011829323800546777, 0.798902627, 0.04333345, 0.358434298, 0.012610262126794884, 0.741299424, 0.0466667, 0.386166944, 0.013391200453042984, 0.6837, 0.05, 0.4139]
veloRekoLUT.ColorSpace = 'RGB'
veloRekoLUT.NanColor = [1.0, 0.0, 0.0]
veloRekoLUT.ScalarRangeInitialized = 1.0
veloRekoLUT.VectorComponent = 2
veloRekoLUT.VectorMode = 'Component'

# trace defaults for the display properties.
veloVectorsDisplay.Representation = 'Surface'
veloVectorsDisplay.ColorArrayName = ['POINTS', veloField]
veloVectorsDisplay.LookupTable = veloRekoLUT
veloVectorsDisplay.Opacity = 0.48
veloVectorsDisplay.SelectTCoordArray = 'None'
veloVectorsDisplay.SelectNormalArray = 'None'
veloVectorsDisplay.SelectTangentArray = 'None'
veloVectorsDisplay.OSPRayScaleArray = veloField
veloVectorsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
veloVectorsDisplay.SelectOrientationVectors = veloField
veloVectorsDisplay.ScaleFactor = 0.06399999856948853
veloVectorsDisplay.SelectScaleArray = 'None'
veloVectorsDisplay.GlyphType = 'Arrow'
veloVectorsDisplay.GlyphTableIndexArray = 'None'
veloVectorsDisplay.GaussianRadius = 0.0031999999284744265
veloVectorsDisplay.SetScaleArray = ['POINTS', veloField]
veloVectorsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
veloVectorsDisplay.OpacityArray = ['POINTS', veloField]
veloVectorsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
veloVectorsDisplay.DataAxesGrid = 'GridAxesRepresentation'
veloVectorsDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
veloVectorsDisplay.ScaleTransferFunction.Points = [-0.0006644729874096811, 0.0, 0.5, 0.0, 0.000712149019818753, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
veloVectorsDisplay.OpacityTransferFunction.Points = [-0.0006644729874096811, 0.0, 0.5, 0.0, 0.000712149019818753, 1.0, 0.5, 0.0]

# show data from annotateTime1
annotateTime1Display = Show(annotateTime1, renderView1, 'TextSourceRepresentation')

# show data from streamTracer1
streamTracer1Display = Show(streamTracer1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
streamTracer1Display.Representation = 'Surface'
streamTracer1Display.ColorArrayName = [None, '']
streamTracer1Display.SelectTCoordArray = 'None'
streamTracer1Display.SelectNormalArray = 'None'
streamTracer1Display.SelectTangentArray = 'None'
streamTracer1Display.OSPRayScaleArray = 'AngularVelocity'
streamTracer1Display.OSPRayScaleFunction = 'PiecewiseFunction'
streamTracer1Display.SelectOrientationVectors = 'Normals'
streamTracer1Display.ScaleFactor = 0.06399999856948853
streamTracer1Display.SelectScaleArray = 'AngularVelocity'
streamTracer1Display.GlyphType = 'Arrow'
streamTracer1Display.GlyphTableIndexArray = 'AngularVelocity'
streamTracer1Display.GaussianRadius = 0.0031999999284744265
streamTracer1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
streamTracer1Display.ScaleTransferFunction = 'PiecewiseFunction'
streamTracer1Display.OpacityArray = ['POINTS', 'AngularVelocity']
streamTracer1Display.OpacityTransferFunction = 'PiecewiseFunction'
streamTracer1Display.DataAxesGrid = 'GridAxesRepresentation'
streamTracer1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
streamTracer1Display.ScaleTransferFunction.Points = [-0.35046177557911085, 0.0, 0.5, 0.0, 0.24977773658914393, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
streamTracer1Display.OpacityTransferFunction.Points = [-0.35046177557911085, 0.0, 0.5, 0.0, 0.24977773658914393, 1.0, 0.5, 0.0]

# show data from tube1
tube1Display = Show(tube1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tube1Display.Representation = 'Surface'
tube1Display.ColorArrayName = ['POINTS', veloField]
tube1Display.LookupTable = veloRekoLUT
tube1Display.Opacity = 0.35
tube1Display.Specular = 1.0
tube1Display.SelectTCoordArray = 'None'
tube1Display.SelectNormalArray = 'TubeNormals'
tube1Display.SelectTangentArray = 'None'
tube1Display.OSPRayScaleArray = 'AngularVelocity'
tube1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tube1Display.SelectOrientationVectors = 'Normals'
tube1Display.ScaleFactor = 0.06444022655487061
tube1Display.SelectScaleArray = 'AngularVelocity'
tube1Display.GlyphType = 'Arrow'
tube1Display.GlyphTableIndexArray = 'AngularVelocity'
tube1Display.GaussianRadius = 0.0032220113277435306
tube1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
tube1Display.ScaleTransferFunction = 'PiecewiseFunction'
tube1Display.OpacityArray = ['POINTS', 'AngularVelocity']
tube1Display.OpacityTransferFunction = 'PiecewiseFunction'
tube1Display.DataAxesGrid = 'GridAxesRepresentation'
tube1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tube1Display.ScaleTransferFunction.Points = [-0.35046177557911085, 0.0, 0.5, 0.0, 0.24977773658914393, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tube1Display.OpacityTransferFunction.Points = [-0.35046177557911085, 0.0, 0.5, 0.0, 0.24977773658914393, 1.0, 0.5, 0.0]

# show data from streamTracer2
streamTracer2Display = Show(streamTracer2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
streamTracer2Display.Representation = 'Surface'
streamTracer2Display.ColorArrayName = [None, '']
streamTracer2Display.SelectTCoordArray = 'None'
streamTracer2Display.SelectNormalArray = 'None'
streamTracer2Display.SelectTangentArray = 'None'
streamTracer2Display.OSPRayScaleArray = 'AngularVelocity'
streamTracer2Display.OSPRayScaleFunction = 'PiecewiseFunction'
streamTracer2Display.SelectOrientationVectors = 'Normals'
streamTracer2Display.ScaleFactor = 0.06399999856948853
streamTracer2Display.SelectScaleArray = 'AngularVelocity'
streamTracer2Display.GlyphType = 'Arrow'
streamTracer2Display.GlyphTableIndexArray = 'AngularVelocity'
streamTracer2Display.GaussianRadius = 0.0031999999284744265
streamTracer2Display.SetScaleArray = ['POINTS', 'AngularVelocity']
streamTracer2Display.ScaleTransferFunction = 'PiecewiseFunction'
streamTracer2Display.OpacityArray = ['POINTS', 'AngularVelocity']
streamTracer2Display.OpacityTransferFunction = 'PiecewiseFunction'
streamTracer2Display.DataAxesGrid = 'GridAxesRepresentation'
streamTracer2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
streamTracer2Display.ScaleTransferFunction.Points = [-0.1591157858314646, 0.0, 0.5, 0.0, 0.5141507281865857, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
streamTracer2Display.OpacityTransferFunction.Points = [-0.1591157858314646, 0.0, 0.5, 0.0, 0.5141507281865857, 1.0, 0.5, 0.0]

# show data from tube2
tube2Display = Show(tube2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tube2Display.Representation = 'Surface'
tube2Display.ColorArrayName = ['POINTS', veloField]
tube2Display.LookupTable = veloRekoLUT
tube2Display.Opacity = 0.35
tube2Display.Specular = 1.0
tube2Display.SelectTCoordArray = 'None'
tube2Display.SelectNormalArray = 'TubeNormals'
tube2Display.SelectTangentArray = 'None'
tube2Display.OSPRayScaleArray = 'AngularVelocity'
tube2Display.OSPRayScaleFunction = 'PiecewiseFunction'
tube2Display.SelectOrientationVectors = 'Normals'
tube2Display.ScaleFactor = 0.06443229019641876
tube2Display.SelectScaleArray = 'AngularVelocity'
tube2Display.GlyphType = 'Arrow'
tube2Display.GlyphTableIndexArray = 'AngularVelocity'
tube2Display.GaussianRadius = 0.0032216145098209383
tube2Display.SetScaleArray = ['POINTS', 'AngularVelocity']
tube2Display.ScaleTransferFunction = 'PiecewiseFunction'
tube2Display.OpacityArray = ['POINTS', 'AngularVelocity']
tube2Display.OpacityTransferFunction = 'PiecewiseFunction'
tube2Display.DataAxesGrid = 'GridAxesRepresentation'
tube2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tube2Display.ScaleTransferFunction.Points = [-0.15128160477060348, 0.0, 0.5, 0.0, 0.499516795139901, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tube2Display.OpacityTransferFunction.Points = [-0.15128160477060348, 0.0, 0.5, 0.0, 0.499516795139901, 1.0, 0.5, 0.0]

# show data from streamTracer3
streamTracer3Display = Show(streamTracer3, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
streamTracer3Display.Representation = 'Surface'
streamTracer3Display.ColorArrayName = [None, '']
streamTracer3Display.SelectTCoordArray = 'None'
streamTracer3Display.SelectNormalArray = 'None'
streamTracer3Display.SelectTangentArray = 'None'
streamTracer3Display.OSPRayScaleArray = 'AngularVelocity'
streamTracer3Display.OSPRayScaleFunction = 'PiecewiseFunction'
streamTracer3Display.SelectOrientationVectors = 'Normals'
streamTracer3Display.ScaleFactor = 0.06399999856948853
streamTracer3Display.SelectScaleArray = 'AngularVelocity'
streamTracer3Display.GlyphType = 'Arrow'
streamTracer3Display.GlyphTableIndexArray = 'AngularVelocity'
streamTracer3Display.GaussianRadius = 0.0031999999284744265
streamTracer3Display.SetScaleArray = ['POINTS', 'AngularVelocity']
streamTracer3Display.ScaleTransferFunction = 'PiecewiseFunction'
streamTracer3Display.OpacityArray = ['POINTS', 'AngularVelocity']
streamTracer3Display.OpacityTransferFunction = 'PiecewiseFunction'
streamTracer3Display.DataAxesGrid = 'GridAxesRepresentation'
streamTracer3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
streamTracer3Display.ScaleTransferFunction.Points = [-0.2822827511038092, 0.0, 0.5, 0.0, 0.21910108921318583, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
streamTracer3Display.OpacityTransferFunction.Points = [-0.2822827511038092, 0.0, 0.5, 0.0, 0.21910108921318583, 1.0, 0.5, 0.0]

# show data from tube3
tube3Display = Show(tube3, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tube3Display.Representation = 'Surface'
tube3Display.ColorArrayName = ['POINTS', veloField]
tube3Display.LookupTable = veloRekoLUT
tube3Display.Opacity = 0.35
tube3Display.Specular = 1.0
tube3Display.SelectTCoordArray = 'None'
tube3Display.SelectNormalArray = 'TubeNormals'
tube3Display.SelectTangentArray = 'None'
tube3Display.OSPRayScaleArray = 'AngularVelocity'
tube3Display.OSPRayScaleFunction = 'PiecewiseFunction'
tube3Display.SelectOrientationVectors = 'Normals'
tube3Display.ScaleFactor = 0.06344650387763977
tube3Display.SelectScaleArray = 'AngularVelocity'
tube3Display.GlyphType = 'Arrow'
tube3Display.GlyphTableIndexArray = 'AngularVelocity'
tube3Display.GaussianRadius = 0.0031723251938819887
tube3Display.SetScaleArray = ['POINTS', 'AngularVelocity']
tube3Display.ScaleTransferFunction = 'PiecewiseFunction'
tube3Display.OpacityArray = ['POINTS', 'AngularVelocity']
tube3Display.OpacityTransferFunction = 'PiecewiseFunction'
tube3Display.DataAxesGrid = 'GridAxesRepresentation'
tube3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tube3Display.ScaleTransferFunction.Points = [-0.2822827511038092, 0.0, 0.5, 0.0, 0.21910108921318583, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tube3Display.OpacityTransferFunction.Points = [-0.2822827511038092, 0.0, 0.5, 0.0, 0.21910108921318583, 1.0, 0.5, 0.0]

# show data from streamTracer4
streamTracer4Display = Show(streamTracer4, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
streamTracer4Display.Representation = 'Surface'
streamTracer4Display.ColorArrayName = [None, '']
streamTracer4Display.SelectTCoordArray = 'None'
streamTracer4Display.SelectNormalArray = 'None'
streamTracer4Display.SelectTangentArray = 'None'
streamTracer4Display.OSPRayScaleArray = 'AngularVelocity'
streamTracer4Display.OSPRayScaleFunction = 'PiecewiseFunction'
streamTracer4Display.SelectOrientationVectors = 'Normals'
streamTracer4Display.ScaleFactor = 0.06399999856948853
streamTracer4Display.SelectScaleArray = 'AngularVelocity'
streamTracer4Display.GlyphType = 'Arrow'
streamTracer4Display.GlyphTableIndexArray = 'AngularVelocity'
streamTracer4Display.GaussianRadius = 0.0031999999284744265
streamTracer4Display.SetScaleArray = ['POINTS', 'AngularVelocity']
streamTracer4Display.ScaleTransferFunction = 'PiecewiseFunction'
streamTracer4Display.OpacityArray = ['POINTS', 'AngularVelocity']
streamTracer4Display.OpacityTransferFunction = 'PiecewiseFunction'
streamTracer4Display.DataAxesGrid = 'GridAxesRepresentation'
streamTracer4Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
streamTracer4Display.ScaleTransferFunction.Points = [-0.2422952550233433, 0.0, 0.5, 0.0, 0.2822477094384668, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
streamTracer4Display.OpacityTransferFunction.Points = [-0.2422952550233433, 0.0, 0.5, 0.0, 0.2822477094384668, 1.0, 0.5, 0.0]

# show data from tube4
tube4Display = Show(tube4, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tube4Display.Representation = 'Surface'
tube4Display.ColorArrayName = ['POINTS', veloField]
tube4Display.LookupTable = veloRekoLUT
tube4Display.Opacity = 0.35
tube4Display.Specular = 1.0
tube4Display.SelectTCoordArray = 'None'
tube4Display.SelectNormalArray = 'TubeNormals'
tube4Display.SelectTangentArray = 'None'
tube4Display.OSPRayScaleArray = 'AngularVelocity'
tube4Display.OSPRayScaleFunction = 'PiecewiseFunction'
tube4Display.SelectOrientationVectors = 'Normals'
tube4Display.ScaleFactor = 0.06493109464645386
tube4Display.SelectScaleArray = 'AngularVelocity'
tube4Display.GlyphType = 'Arrow'
tube4Display.GlyphTableIndexArray = 'AngularVelocity'
tube4Display.GaussianRadius = 0.003246554732322693
tube4Display.SetScaleArray = ['POINTS', 'AngularVelocity']
tube4Display.ScaleTransferFunction = 'PiecewiseFunction'
tube4Display.OpacityArray = ['POINTS', 'AngularVelocity']
tube4Display.OpacityTransferFunction = 'PiecewiseFunction'
tube4Display.DataAxesGrid = 'GridAxesRepresentation'
tube4Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tube4Display.ScaleTransferFunction.Points = [-0.2422952550233433, 0.0, 0.5, 0.0, 0.2822477094384668, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tube4Display.OpacityTransferFunction.Points = [-0.2422952550233433, 0.0, 0.5, 0.0, 0.2822477094384668, 1.0, 0.5, 0.0]

# show data from vortexCores1
vortexCores1Display = Show(vortexCores1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
vortexCores1Display.Representation = 'Surface'
vortexCores1Display.ColorArrayName = [None, '']
vortexCores1Display.SelectTCoordArray = 'None'
vortexCores1Display.SelectNormalArray = 'None'
vortexCores1Display.SelectTangentArray = 'None'
vortexCores1Display.OSPRayScaleArray = veloField
vortexCores1Display.OSPRayScaleFunction = 'PiecewiseFunction'
vortexCores1Display.SelectOrientationVectors = 'acceleration'
vortexCores1Display.ScaleFactor = 0.05866659879684449
vortexCores1Display.SelectScaleArray = 'None'
vortexCores1Display.GlyphType = 'Arrow'
vortexCores1Display.GlyphTableIndexArray = 'None'
vortexCores1Display.GaussianRadius = 0.002933329939842224
vortexCores1Display.SetScaleArray = ['POINTS', veloField]
vortexCores1Display.ScaleTransferFunction = 'PiecewiseFunction'
vortexCores1Display.OpacityArray = ['POINTS', veloField]
vortexCores1Display.OpacityTransferFunction = 'PiecewiseFunction'
vortexCores1Display.DataAxesGrid = 'GridAxesRepresentation'
vortexCores1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
vortexCores1Display.ScaleTransferFunction.Points = [-0.0020564887672662735, 0.0, 0.5, 0.0, 0.0008443798869848251, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
vortexCores1Display.OpacityTransferFunction.Points = [-0.0020564887672662735, 0.0, 0.5, 0.0, 0.0008443798869848251, 1.0, 0.5, 0.0]

# show data from tube5
tube5Display = Show(tube5, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tube5Display.Representation = 'Surface'
tube5Display.AmbientColor = [0.0, 0.0, 0.0]
tube5Display.ColorArrayName = [None, '']
tube5Display.DiffuseColor = [0.0, 0.0, 0.0]
tube5Display.SelectTCoordArray = 'None'
tube5Display.SelectNormalArray = 'TubeNormals'
tube5Display.SelectTangentArray = 'None'
tube5Display.OSPRayScaleArray = 'TubeNormals'
tube5Display.OSPRayScaleFunction = 'PiecewiseFunction'
tube5Display.SelectOrientationVectors = 'acceleration'
tube5Display.ScaleFactor = 0.05907695591449738
tube5Display.SelectScaleArray = 'None'
tube5Display.GlyphType = 'Arrow'
tube5Display.GlyphTableIndexArray = 'None'
tube5Display.GaussianRadius = 0.0029538477957248687
tube5Display.SetScaleArray = ['POINTS', 'TubeNormals']
tube5Display.ScaleTransferFunction = 'PiecewiseFunction'
tube5Display.OpacityArray = ['POINTS', 'TubeNormals']
tube5Display.OpacityTransferFunction = 'PiecewiseFunction'
tube5Display.DataAxesGrid = 'GridAxesRepresentation'
tube5Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tube5Display.ScaleTransferFunction.Points = [-0.9939237236976624, 0.0, 0.5, 0.0, 0.9939237236976624, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tube5Display.OpacityTransferFunction.Points = [-0.9939237236976624, 0.0, 0.5, 0.0, 0.9939237236976624, 1.0, 0.5, 0.0]

# show data from bfieldVektor
bfieldVektorDisplay = Show(bfieldVektor, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
bfieldVektorDisplay.Representation = 'Surface'
bfieldVektorDisplay.ColorArrayName = ['POINTS', 'bfield_sensor']
bfieldVektorDisplay.LookupTable = bfield_sensorLUT
bfieldVektorDisplay.SelectTCoordArray = 'None'
bfieldVektorDisplay.SelectNormalArray = 'None'
bfieldVektorDisplay.SelectTangentArray = 'None'
bfieldVektorDisplay.OSPRayScaleArray = 'bfield_sensor'
bfieldVektorDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
bfieldVektorDisplay.SelectOrientationVectors = 'sensorNormals'
bfieldVektorDisplay.ScaleFactor = 0.5137511312961579
bfieldVektorDisplay.SelectScaleArray = 'bfield_sensor'
bfieldVektorDisplay.GlyphType = 'Arrow'
bfieldVektorDisplay.GlyphTableIndexArray = 'bfield_sensor'
bfieldVektorDisplay.GaussianRadius = 0.02568755656480789
bfieldVektorDisplay.SetScaleArray = ['POINTS', 'bfield_sensor']
bfieldVektorDisplay.ScaleTransferFunction = 'PiecewiseFunction'
bfieldVektorDisplay.OpacityArray = ['POINTS', 'bfield_sensor']
bfieldVektorDisplay.OpacityTransferFunction = 'PiecewiseFunction'
bfieldVektorDisplay.DataAxesGrid = 'GridAxesRepresentation'
bfieldVektorDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
bfieldVektorDisplay.ScaleTransferFunction.Points = [-57.95941162109375, 0.0, 0.5, 0.0, 83.021728515625, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
bfieldVektorDisplay.OpacityTransferFunction.Points = [-57.95941162109375, 0.0, 0.5, 0.0, 83.021728515625, 1.0, 0.5, 0.0]

# show data from sensorBox
sensorBoxDisplay = Show(sensorBox, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
sensorBoxDisplay.Representation = 'Surface'
sensorBoxDisplay.ColorArrayName = ['POINTS', 'bfield_sensor']
sensorBoxDisplay.LookupTable = bfield_sensorLUT
sensorBoxDisplay.SelectTCoordArray = 'None'
sensorBoxDisplay.SelectNormalArray = 'Normals'
sensorBoxDisplay.SelectTangentArray = 'None'
sensorBoxDisplay.OSPRayScaleArray = 'bfield_sensor'
sensorBoxDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
sensorBoxDisplay.SelectOrientationVectors = 'sensorNormals'
sensorBoxDisplay.ScaleFactor = 0.05
sensorBoxDisplay.SelectScaleArray = 'bfield_sensor'
sensorBoxDisplay.GlyphType = 'Arrow'
sensorBoxDisplay.GlyphTableIndexArray = 'bfield_sensor'
sensorBoxDisplay.GaussianRadius = 0.0025
sensorBoxDisplay.SetScaleArray = ['POINTS', 'bfield_sensor']
sensorBoxDisplay.ScaleTransferFunction = 'PiecewiseFunction'
sensorBoxDisplay.OpacityArray = ['POINTS', 'bfield_sensor']
sensorBoxDisplay.OpacityTransferFunction = 'PiecewiseFunction'
sensorBoxDisplay.DataAxesGrid = 'GridAxesRepresentation'
sensorBoxDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
sensorBoxDisplay.ScaleTransferFunction.Points = [-57.95941162109375, 0.0, 0.5, 0.0, 83.021728515625, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sensorBoxDisplay.OpacityTransferFunction.Points = [-57.95941162109375, 0.0, 0.5, 0.0, 83.021728515625, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for bfield_sensorLUT in view renderView1
bfield_sensorLUTColorBar = GetScalarBar(bfield_sensorLUT, renderView1)
bfield_sensorLUTColorBar.Title = 'bfield in nT'
bfield_sensorLUTColorBar.ComponentTitle = ''

# set color bar visibility
bfield_sensorLUTColorBar.Visibility = 0

# get color legend/bar for temperatureLUT in view renderView1
temperatureLUTColorBar = GetScalarBar(temperatureLUT, renderView1)
temperatureLUTColorBar.WindowLocation = 'Any Location'
temperatureLUTColorBar.Position = [0.013636363636363606, 0.5851197982345523]
temperatureLUTColorBar.Title = 'temperature in Â°C'
temperatureLUTColorBar.ComponentTitle = ''
temperatureLUTColorBar.ScalarBarLength = 0.2719924337957126

# set color bar visibility
temperatureLUTColorBar.Visibility = 1

# get color legend/bar for veloRekoLUT in view renderView1
veloRekoLUTColorBar = GetScalarBar(veloRekoLUT, renderView1)
veloRekoLUTColorBar.WindowLocation = 'Any Location'
veloRekoLUTColorBar.Position = [0.019480519480519445, 0.21311475409836061]
veloRekoLUTColorBar.Title = 'VeloReko in m/s'
veloRekoLUTColorBar.ComponentTitle = 'Z'
veloRekoLUTColorBar.ScalarBarLength = 0.2682093316519545

# set color bar visibility
veloRekoLUTColorBar.Visibility = 1

# hide data in view
Hide(bfield_sensor_velo_, renderView1)

# show color legend
temp___Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(temp___, renderView1)

# show color legend
tempBoxesDisplay.SetScalarBarVisibility(renderView1, True)

# show color legend
veloVectorsDisplay.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(streamTracer1, renderView1)

# show color legend
tube1Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(tube1, renderView1)

# hide data in view
Hide(streamTracer2, renderView1)

# show color legend
tube2Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(tube2, renderView1)

# hide data in view
Hide(streamTracer3, renderView1)

# show color legend
tube3Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(tube3, renderView1)

# hide data in view
Hide(streamTracer4, renderView1)

# show color legend
tube4Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(tube4, renderView1)

# hide data in view
Hide(vortexCores1, renderView1)

# hide data in view
Hide(bfieldVektor, renderView1)

# hide data in view
Hide(sensorBox, renderView1)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'bfield_sensor'
bfield_sensorPWF = GetOpacityTransferFunction('bfield_sensor')
bfield_sensorPWF.Points = [-57.95941162109375, 0.0, 0.5, 0.0, 83.021728515625, 1.0, 0.5, 0.0]
bfield_sensorPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for veloField
veloRekoPWF = GetOpacityTransferFunction(veloField)
veloRekoPWF.Points = [-0.01940820924937725, 0.0, 0.5, 0.0, 0.013391200453042984, 1.0, 0.5, 0.0]
veloRekoPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'temperature'
temperaturePWF = GetOpacityTransferFunction('temperature')
temperaturePWF.Points = [22.27692985534668, 0.0, 0.5, 0.0, 22.503889083862305, 1.0, 0.5, 0.0]
temperaturePWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# restore active source
SetActiveSource(tube5)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
   pass
   
