#import various packages needed to generate decision tree 
import numpy,Slope_aspect,city_block_dist,time,rpy,pickle
def DecisionTree():
  #importing rpart from rpy
  rpy.r.library("rpart")

  time_a = time.time()
  print "Reading Elevation Data from ascii file"
  Elev_arr        = numpy.loadtxt("training_data/new_elev.asc",unpack=True)
  time_b = time.time()

  print "Reading Landcover Data from ascii file" , (time_b - time_a)
  Landcover      = numpy.loadtxt("training_data/new_land.asc",unpack=True)
  time_b = time.time()

  print "Reading River Data from ascii file" , (time_b - time_a)
  River          = numpy.loadtxt("training_data/new_rivers.asc",unpack=True)
  time_b = time.time()

  print "Computing City block distance from River data", (time_b - time_a)
  River_dist_arr = city_block_dist.CityBlock(River)
  time_b = time.time()

  print "Computing Slope and Aspect from Elevation data", (time_b - time_a)
  (Slope_arr,Aspect_arr) = Slope_aspect.Slope_aspect(Elev_arr)
  time_b = time.time()

  (x_len,y_len) = Elev_arr.shape
  no_of_veg_class = 10 #no of vegetation class in Landcover matrix
  print "Generating Lists of differnt Landcover classes", (time_b - time_a)
 
  #Create list of lists to hold pixels of each landcover class - no of list in 
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

  time_b = time.time()
  print "Sampling Data for decision tree" , (time_b - time_a)

  #normalizing elevation data
  #minimum_elev = numpy.min(Elev_arr)
  #factor = numpy.max(Elev_arr) - minimum_elev
  #Elev_arr = (Elev_arr[:,:]-minimum_elev)*100/factor

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
  rpy.r.png('Output/DecisionTree.png')
  rpy.r.plot(fit)
  rpy.r.text(fit)
  rpy.r.dev_off()
