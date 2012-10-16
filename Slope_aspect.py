import numpy,scipy

def Slope_aspect(DEM_arr):
  """
  Return the Slope and Aspect given the Digital Elevation Model
  """
  (x_len,y_len) = DEM_arr.shape
  
  # Assigning different code to different direction
  # based on some logic described below
  #
  #   N-W (-1, 1)        North(0,1)      N-E(1,1)
  #                        |
  #                        |
  #   West(-1, 0)__________|___________  East(1,0)
  #                        |
  #                        |
  #   S-W (-1,-1)        South(0,-1)      S-E(1,-1) 

  Flat       = (0,0)
  North      = (0,1)
  North_East = (1,1)
  East       = (1,0)
  South_East = (1,-1)
  South      = (0,-1)
  South_West = (-1,-1)
  West       = (-1,0)
  North_West = (-1,1)

  #Allocating array to hold Slope and Aspect values
  Aspect_arr = numpy.zeros((x_len,y_len,2),dtype ="int" )
  slope_arr = numpy.zeros((x_len,y_len),dtype ="uint16" )

  for i in range(1,x_len-1):
    for j in range(1,y_len-1): #Boundary cells have been ignored
      
      #Formula for calculating Slope and Aspect is same as in ArcGIS.so, have a 
      # look on the documentation for any doubt.
      rate_change_x =((DEM_arr[i-1][j+1]+2*DEM_arr[i][j+1]+DEM_arr[i+1][j+1])\
                     -(DEM_arr[i-1][j-1]+2*DEM_arr[i][j-1]+DEM_arr[i+1][j-1]))/8
      rate_change_y =((DEM_arr[i+1][j-1]+2*DEM_arr[i+1][j]+DEM_arr[i+1][j+1])\
                     -(DEM_arr[i-1][j-1]+2*DEM_arr[i-1][j]+DEM_arr[i-1][j+1]))/8
      
      #_______________________Aspect Calculation__________________________________ 
      aspect = 57.29578*numpy.arctan2(rate_change_y,-1*rate_change_x)
      if aspect == 0.0:
        cell = 361 # to handle pixel in flat area
      elif aspect < 0:
        cell = 90.0 - aspect
      elif aspect > 90.0:
        cell = 360.0 - aspect + 90.0 
      else:
        cell = 90.0 - aspect

      if cell == 361:
        Aspect_arr[i,j] = Flat       # Flat land
      elif  (0.0 <= cell < 22.) or (337.5 <= cell < 360.0):
        Aspect_arr[i,j] = North      # North
      elif ( 22.5 <= cell < 67.5 ):
        Aspect_arr[i,j] = North_East # North-East
      elif ( 67.5 <= cell < 112.5 ):
        Aspect_arr[i,j] = East       # East
      elif ( 112.5 <= cell < 157.5 ):
        Aspect_arr[i,j] = South_East # South-East
      elif ( 157.5 <= cell < 202.5 ):
        Aspect_arr[i,j] = South      # South
      elif ( 202.5 <= cell < 247.5 ):
        Aspect_arr[i,j] = South_West # South-West
      elif ( 247.5 <= cell < 292.5 ):
        Aspect_arr[i,j] = West       # West
      elif ( 292.5 <= cell < 337.5 ):
        Aspect_arr[i,j] = North_West # North-West

      #_____________________Slope Calculation_____________________________________    
      rise_run = numpy.sqrt((rate_change_x*rate_change_x + rate_change_y*rate_change_y)/900)
      slope_degree = 57.29*numpy.arctan(rise_run)
      slope_arr[i][j] = slope_degree

  return (slope_arr,Aspect_arr)
