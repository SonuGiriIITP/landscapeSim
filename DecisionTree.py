import numpy
import Slope_aspect
import city_block_dist
import rpy
import pickle

def DecisionTree(output_dir, elev_filename, landcover_filename, river_filename):
  """
  This module generate decision tree used to allocate landcover classes.
  It imports rpart library from rpy package.
  Reads the training data, creates a sample data and use rpart libray to build decision tree. 
  """
  rpy.r.library("rpart") # rpart library used for creating Decision tree
  #Read Elevation Data from ascii file
  file_name = "training_data/%s" % (elev_filename)
  Elev_arr = numpy.loadtxt(file_name, unpack=True)
  #Read Landcover Data from ascii file
  file_name = "training_data/%s" % (landcover_filename)
  Landcover = numpy.loadtxt(file_name, unpack=True)
  #Read River Data from ascii file
  file_name = "training_data/%s" % (river_filename)
  River     = numpy.loadtxt(file_name, unpack=True)
  #Compute City block distance from River data
  River_dist_arr = city_block_dist.CityBlock(River)
  #Compute Slope and Aspect from Elevation data
  (Slope_arr,Aspect_arr) = Slope_aspect.Slope_aspect(Elev_arr)

  (x_len,y_len) = Elev_arr.shape
  no_of_veg_class = 10 #no of vegetation class in Landcover matrix
  #Generating Lists for differnt Landcover classes
  # Create list of lists to hold pixels of each landcover class - no of list in 
  # list L is equal to no_of_veg_class  
  L = []
  for i in range(0,no_of_veg_class):
    L.append([]) 

  #Now append the pixel co-ordinates into respective list of lists
  for i in range(1,x_len-1):   # Ignoring boundary cells 
    for j in range(1,y_len-1): # because we don't have slope and aspect for them
      #nodata values already gets handled since we are ignoring it
      if Landcover[i][j] == 0:
        L[0].append( ( i,j ) )
      elif Landcover[i][j] == 1:
        L[1].append( ( i,j ) )
      elif Landcover[i][j] == 2:
        L[2].append( ( i,j ) )
      elif Landcover[i][j] == 3:
        L[3].append( ( i,j ) )
      elif Landcover[i][j] == 4:
        L[4].append( ( i,j ) )
      elif Landcover[i][j] == 5:
        L[5].append( ( i,j ) )
      elif Landcover[i][j] == 6:
        L[6].append( ( i,j ) )
      elif Landcover[i][j] == 7:
        L[7].append( ( i,j ) )
      elif Landcover[i][j] == 8:
        L[8].append( ( i,j ) )
      elif Landcover[i][j] == 9:
        L[9].append( ( i,j ) )

  #Sample Data for decision tree
  #normalizing elevation data
  minimum_elev = numpy.min(Elev_arr)
  factor = numpy.max(Elev_arr) - minimum_elev
  Elev_arr = (Elev_arr[:,:]-minimum_elev)*100/factor

  #Create various list to hold sample training data
  Elevation = []
  Slope = []
  RiverDistance = []
  Aspect_x = []
  Aspect_y = []
  Class = []

  #Now sampling the data
  for i in range(0,no_of_veg_class):  
    if len(L[i]) < 500:
      limit = len(L[i])
    else:
      limit = 500
    for j in range(0,limit):
      Elevation.append( int(Elev_arr[ L[i][j][0] ][ L[i][j][1] ]))
      Slope.append(int(Slope_arr[ L[i][j][0] ][ L[i][j][1] ]))
      RiverDistance.append(int(River_dist_arr[ L[i][j][0] ][ L[i][j][1] ]))
      Aspect_x.append(int(Aspect_arr[ L[i][j][0] ][ L[i][j][1] ][0]))
      Aspect_y.append(int(Aspect_arr[ L[i][j][0] ][ L[i][j][1] ][1]))
      Class.append(i)

  #create dictionary of sample data which will be needed to generate decision tree 
  traing_data = {'Elevation':Elevation,'Slope':Slope,'RiverDistance':RiverDistance,'Aspect_x':Aspect_x,'Aspect_y':Aspect_y,'Class':Class}

  #write dictionary into pickle file for further use(reusability) 
  output =  open('decision_tree.pkl','wb')
  pickle.dump(traing_data,output)
  output.close()

  rpy.set_default_mode(rpy.NO_CONVERSION)
  print "Creating Decision tree"
  #Using rpart create the decision tree
  fit = rpy.r.rpart(formula='Class ~ Elevation + RiverDistance + Slope + Aspect_x + Aspect_y',data=traing_data,method="class")

  #output a png image of the decision tree
  file_name = "%s/DecisionTree.png" % (output_dir)
  rpy.r.png(file_name)
  rpy.r.plot(fit)
  rpy.r.text(fit)
  rpy.r.dev_off()
