import numpy
import Slope_aspect
import dist
import time
import rpy
import pickle
import pylab
from scipy import ndimage

def VegetationClassify(Elev_arr, River_arr): 

  rpy.r.library("rpart")
  # Read the dictionary from the pickle file
  pkl_file = open('decision_tree.pkl','rb')
  rpy.set_default_mode(rpy.NO_CONVERSION)
  traing_data = pickle.load(pkl_file)
  pkl_file.close()

  # Create Decision tree for predicting landcover class
  # create the decision tree using rpart 
  fit = rpy.r.rpart(formula='Class ~ Elevation + RiverDistance + Slope \
      + Aspect_x + Aspect_y',data = traing_data, method = "class")

  # calculate River distance using River_arr
  River_dist_arr = dist.CityBlock(River_arr)  
  # claculate slope and aspect
  (Slope_arr, Aspect_arr) = Slope_aspect.Slope_aspect(Elev_arr)

  (x_len, y_len) = Elev_arr.shape
  # Alloctae vegetation array for holding predicted landcover values
  Veg_arr = numpy.zeros((x_len, y_len), dtype = "uint8")

  # Normalize the elevation data
  minimum_elev = numpy.min(Elev_arr)
  factor = numpy.max(Elev_arr) - minimum_elev
  Elev_arr = (Elev_arr[:,:] - minimum_elev)*100/factor

  # Create various list to hold test data
  Elevation = []
  Slope = []
  RiverDistance = []
  Aspect_x = []
  Aspect_y = []

  # Append the data into respective list
  for i in range(0,x_len):
    for j in range(0,y_len):
      Elevation.append(int(Elev_arr[i][j]))
      Slope.append(int(Slope_arr[i][j]))
      RiverDistance.append(int(River_dist_arr[i][j]))
      Aspect_x.append(int(Aspect_arr[i][j][0]))
      Aspect_y.append(int(Aspect_arr[i][j][1]))

  # Create dictionary so as to apply R's predict command on it 
  Test_data ={'Elevation':Elevation ,'Slope':Slope ,'RiverDistance':RiverDistance,\
             'Aspect_x':Aspect_x,'Aspect_y':Aspect_y}

  rpy.set_default_mode(rpy.BASIC_CONVERSION)
  # values contain probability values of the predicted landcover classes
  values = rpy.r.predict(fit,newdata=Test_data,method="class")
  for i in range(0,x_len):
    for j in range(0,y_len):
      # Get the class having max probability for each test data point
      a = ndimage.maximum_position(values[i*x_len + j])
      Veg_arr[i,j] = (a[0]*25) # Assign them some value to facilitate visualization
  return Veg_arr
