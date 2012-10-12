import Hydro_Network
import DecisionTree
import VegetationClassify 
import Geometry
import time

raw_input("Press Enter to continue...else exit usind CTRL + Z")
print "Please wait...while the program is executing"
print "Running Digital Elevation Model and River Network module"
time1 = time.time()

Hydro_Network.RiverNetwork(H1 = 0.9, H2 = 0.85 , H3 = 0.70, Elev_range = [0,1309])
time2 = time.time()
print "Time taken to generate DEM , River Network, Catchment Matrix is" , (time2 -time1),"seconds"

response = raw_input("Press y/Y if you want to use training data to build decision tree ,press any other key to skip")
if (response == 'y') or (response == 'Y'):
  DecisionTree.DecisionTree()
  time3 = time.time()
  print "Time taken to generate decision tree is " , (time3 - time2) ,"seconds"

time3 = time.time()
VegetationClassify.VegetationClassify()
time4 = time.time()
print "Time taken to assign landcover is " , (time4 - time3),"seconds"

Geometry.GeometricFeature(min_area = 40,max_area = 400,aspect_ratio = 1.8,agri_area_limit = 0.3)
time5 = time.time()
print "Time taken to generate Geometric Features is " ,(time5 - time4) ,"seconds"
