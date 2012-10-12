import numpy,scipy,math,pylab,MapGeneration_pure_python,Flow_accum
from scipy import ndimage
import city_block_dist_erosion

def pixel_exist(p,q,x_len,y_len):
  """Check weather the pixel lie within the boundary of the grid"""
  return ( (p >= 1) and (q >= 1) and ( p <= x_len - 2) and (q <= x_len-2) )

def RiverNetwork(H1 = 0.9, H2 = 0.85 , H3 = 0.70, Elev_range = [0,1309]):

  print "Generating Digital Elevation Maps using FM2D algorithm"

  #Generate first DEM with gradient = 1 (i.e. TRUE) and high H value 
  DEM_arr1 = MapGeneration_pure_python.midPointFm2d(max_level = 9, sigma = 1, H = H1, addition = True,\
                        wrap = False, gradient = 1,seed = 0, normalise=True,bounds=Elev_range)
  pylab.imsave("Output/DigitalElevationModel1",DEM_arr1)

  #Generate second DEM with gradient = 0 (i.e. FLASE) and medium H value 
  DEM_arr2 = MapGeneration_pure_python.midPointFm2d(max_level = 9, sigma = 1, H = H2, addition = True,\
                        wrap = False, gradient = 0,seed = 65, normalise = True,bounds=Elev_range)
  pylab.imsave("Output/DigitalElevationModel2",DEM_arr2)

  #Generate third DEM with gradient = 0 (i.e. FLASE) and medium H value 
  DEM_arr3 = MapGeneration_pure_python.midPointFm2d(max_level = 9, sigma = 1, H = H3, addition = True,\
                        wrap = False, gradient = 0,seed = 6, normalise = True,bounds=Elev_range)
  pylab.imsave("Output/DigitalElevationModel3",DEM_arr3)

  DEM_arr = DEM_arr1
  (x_len,y_len) = DEM_arr.shape
  #Get the co-ordinates having highest elev , required for catchment extraction
  (max_x, max_y) = ndimage.maximum_position(DEM_arr)

  print "Iteratively removing sink using 3x3 window"
  for p in range(0,6):
    for i in range(1,x_len-1):
      for j in range(1,y_len-1):
        #Remove pits by 3 x 3 window
        A = min(DEM_arr[i-1][j-1],DEM_arr[i-1][j],DEM_arr[i-1][j+1],\
                DEM_arr[i][j-1],DEM_arr[i][j+1],DEM_arr[i+1][j-1],
                DEM_arr[i+1][j],DEM_arr[i+1][j+1])
        if DEM_arr[i][j] < A:
          DEM_arr[i][j] = A + 1

  #Initialize various arrays to hold flow direction, Flow accumulation,catchment info etc
  Flow_arr = numpy.zeros((x_len,y_len) , dtype = "uint8" )
  # River_arr will hold the River_accumulation matrix
  River_arr = numpy.ones((x_len,y_len) , dtype = "int" )
  # Catchment_boundary_arr will hold the Catchment boundaries
  Catchment_boundary_arr = numpy.zeros((x_len,y_len) , dtype = "uint8" )
  # Found_arr will hold the catchment with different labels
  Found_arr = numpy.zeros((x_len,y_len),dtype = "uint8")
  # Pour_point_arr keeps track of pour point of a catchment on the map
  Pour_point_arr = numpy.zeros((x_len,y_len),dtype = "uint8")
  Pour_point_list = [] #keep track of Pour_point in a list

  pit_list = [] #contains all the pit in DEM
  print "Assigning Flow Directions"
  for i in range(1,x_len-1): 
    for j in range(1,y_len-1):
      #Assign Flow direction
      (value,dirn) =max(((DEM_arr[i][j] - DEM_arr[i-1][j-1])/1.41,3),\
                    (DEM_arr[i][j]-DEM_arr[i-1][j],2),((DEM_arr[i][j]-DEM_arr[i-1][j+1])/1.41,1),\
                    (DEM_arr[i][j] - DEM_arr[i][j-1],4),(0,8),(DEM_arr[i][j] - DEM_arr[i][j+1],0),\
                    ((DEM_arr[i][j] - DEM_arr[i+1][j-1])/1.41,5),(DEM_arr[i][j] - DEM_arr[i+1][j],6),\
                    ((DEM_arr[i][j] - DEM_arr[i+1][j+1])/1.41,7))
      Flow_arr[i][j] = dirn
      if dirn == 8:
        # If there is a pit append it to the pit_list
        pit_list.append((i,j))

  label = 0 # will be used to assign labels to differnet catchments

#_____________Catchment Extraction_____________________________________________
  print "Extracting Catchment and filling Depressions"
  while len(pit_list) >= 1:
  #_______________For each and every pit in the DEM do _________________________
    stack = []
    pit = pit_list.pop(0)
    stack.append(pit)
    label = label + 1 #increase the label being assigned to the catchment 
    #_______________________Identify catchment for each and every pit
    catchment_pixels = []
    catchment_pixels.append((DEM_arr[pit[0],pit[1]],pit[0],pit[1]))
    while len(stack) > 0:
      (p,q) = stack.pop(0)
      Found_arr[p][q] = label
      #Pop an element from stack check if its adjacent pixels exist and contribute 
      # its flow to the central pixel(pixel popped) then append it into list, continue
      # this while stack gets empty
      if pixel_exist(p-1,q-1,x_len,y_len):
        if Flow_arr[p-1][q-1] == 7:
          stack.append((p-1,q-1))
          catchment_pixels.append((DEM_arr[p-1,q-1],p-1,q-1))
      if pixel_exist(p-1,q,x_len,y_len):
        if Flow_arr[p-1][q] == 6 :
          catchment_pixels.append((DEM_arr[p-1,q],p-1,q))
          stack.append((p-1,q))
      if pixel_exist(p-1,q+1,x_len,y_len):
        if Flow_arr[p-1][q+1] == 5:
          catchment_pixels.append((DEM_arr[p-1,q+1],p-1,q+1))
          stack.append((p-1,q+1))
      if pixel_exist(p,q-1,x_len,y_len):
        if Flow_arr[p][q-1] == 0 :
          catchment_pixels.append((DEM_arr[p,q-1],p,q-1))
          stack.append((p,q-1))
      if pixel_exist(p,q+1,x_len,y_len):
        if Flow_arr[p][q+1] == 4:
          catchment_pixels.append((DEM_arr[p,q+1],p,q+1))
          stack.append((p,q+1))
      if pixel_exist(p+1,q-1,x_len,y_len):
        if Flow_arr[p+1][q-1] == 1:
          catchment_pixels.append((DEM_arr[p+1,q-1],p+1,q-1))
          stack.append((p+1,q-1))
      if pixel_exist(p+1,q,x_len,y_len):
        if Flow_arr[p+1][q] == 2 :
          catchment_pixels.append((DEM_arr[p+1,q],p+1,q))
          stack.append((p+1,q))
      if pixel_exist(p+1,q+1,x_len,y_len):
        if Flow_arr[p+1][q+1] == 3 :
          catchment_pixels.append((DEM_arr[p+1,q+1],p+1,q+1))
          stack.append((p+1,q+1))
    # Find catchment Outlet
    pour_point = (max_x, max_y)
    flag = 0
    for i in range(0,len(catchment_pixels)):
      (p,q) = ( catchment_pixels[i][1],catchment_pixels[i][2] )
      label = Found_arr[p][q]
      # Catchment Outlet will be the minimum catchment boundary pixel
      if (Found_arr[p-1][q-1] != label or Found_arr[p-1][q] != label or Found_arr[p-1][q+1] != label or 
         Found_arr[p][q-1] != label or Found_arr[p][q+1] != label or Found_arr[p+1][q-1] != label or
         Found_arr[p+1][q] != label or Found_arr[p+1][q+1] != label):# if pixel lie on boundary of catchment
        Catchment_boundary_arr[p][q] = 255
        if DEM_arr[ pour_point[0] ][ pour_point[1] ] > DEM_arr[p][q]:#if height of boundary is less then update pour point
          pour_point = (p,q)
          flag = 1
    if flag == 1:
      Pour_point_list.append((DEM_arr[pour_point],pour_point[0],pour_point[1]))
      Pour_point_arr[pour_point] = 255
      for i in range(0,len(catchment_pixels)):
        if catchment_pixels[i][0] < DEM_arr[pour_point]:
          #fill the depression in the catchment
          DEM_arr[catchment_pixels[i][1],catchment_pixels[i][2]] = DEM_arr[pour_point]

  print "Assignnig flow dirnection again after Depression filling"
  for i in range(1,x_len-1):
    for j in range(1,y_len-1):
    # Again assign Flow direction again after filling the depressions
      (value, dirn ) = max( ((DEM_arr[i][j] - DEM_arr[i-1][j-1])/1.41,3),(DEM_arr[i][j] - DEM_arr[i-1][j],2),((DEM_arr[i][j] - DEM_arr[i-1][j+1])/1.41,1),\
                            (DEM_arr[i][j] - DEM_arr[i][j-1],4),(0,8),(DEM_arr[i][j] - DEM_arr[i][j+1],0),\
                            ((DEM_arr[i][j] - DEM_arr[i+1][j-1])/1.41,5),(DEM_arr[i][j] - DEM_arr[i+1][j],6),((DEM_arr[i][j] - DEM_arr[i+1][j+1])/1.41,7))
      Flow_arr[i][j] = dirn
      if value <= 0:
        Flow_arr[i][j] = 8

  # Calculate flow accumulation by calling Generate_River function
  print "Performing Flow accumulation"
  River_arr  = Flow_accum.Generate_River( Flow_arr,River_arr,DEM_arr)

  Distance_arr = city_block_dist_erosion.CityBlock(River_arr)
  # Create a mask for differnet distances used for DEM erosion
  print "Eroding DEM"
  mask4 = [ Distance_arr <= 15 ]
  mask5 = [ Distance_arr > 3 ]
  mask3 = [Distance_arr == 3]
  mask2 = [Distance_arr == 2]
  mask1 = [Distance_arr == 1]
  mask0 = [Distance_arr == 0]
  max_flow_accum = numpy.max(River_arr)

  for i in range(0,x_len):
    for j in range(0,y_len):
      #Erode the landscape using diffent weighing factor for different distances from 
      #river while combining 3 DEM's 
      if mask0[0][i][j] == True:
        DEM_arr[i][j] = 0.3*DEM_arr[i][j] + 0.45*DEM_arr2[i][j] + 0.19*DEM_arr3[i][j]
      elif mask1[0][i][j] == True:
        DEM_arr[i][j] = 0.3*DEM_arr[i][j] + 0.46*DEM_arr2[i][j] + 0.20*DEM_arr3[i][j]
      elif mask2[0][i][j] == True:
        DEM_arr[i][j] = 0.3*DEM_arr[i][j] + 0.46*DEM_arr2[i][j] + 0.21*DEM_arr3[i][j]
      elif mask3[0][i][j] == True:
        DEM_arr[i][j] = 0.3*DEM_arr[i][j] + 0.46*DEM_arr2[i][j] + 0.23*DEM_arr3[i][j]
      elif mask4[0][i][j] == True and mask5[0][i][j] == True:
        DEM_arr[i][j] = 0.3*DEM_arr[i][j] + 0.47*DEM_arr2[i][j] + 0.23*DEM_arr3[i][j]
      else:     
        DEM_arr[i][j] = 0.3*DEM_arr[i][j] + 0.50*DEM_arr2[i][j] + 0.25*DEM_arr3[i][j]

#Output different statistics for display and further use
  print "printing statistics ...see the Output Folder"
  numpy.save("River.npy",River_arr) 
  numpy.save("DEM.npy",DEM_arr)
  pylab.imsave("Output/River",River_arr)
  pylab.imsave("Output/Catchment",Found_arr)
  pylab.imsave("Output/CatchmentBoundary",Catchment_boundary_arr)
  pylab.imsave("Output/Combined_eroded_DEM",DEM_arr)
  pylab.imsave("Output/RiverDistance",Distance_arr)
