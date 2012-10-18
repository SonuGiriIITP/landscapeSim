import scipy
import numpy
from scipy import *

def CityBlock(River_arr):
  ( x_len , y_len ) = River_arr.shape
  
  #Initialize distance array  
  distance_arr = zeros( (x_len , y_len), dtype= "int")
  infinity = x_len + y_len + 1
  for i in range(0, x_len):
    for j in range(0, y_len):
      # Convert the gray-scale image into binary image
      if River_arr[i][j] == 0: 
        # 0 means absence of river
        distance_arr[i][j] = infinity
      else:
        distance_arr[i][j] = 0

  # Apply distance transform
  for i in range(0, x_len):
    for j in range(0, y_len):
      if i-1 < 0 and j-1 >=0:
        distance_arr[i][j] = min(distance_arr[i][j],1 + distance_arr[i][j-1])
      elif i-1 >= 0 and j-1 < 0:
        distance_arr[i][j] = min(distance_arr[i][j],1 + distance_arr[i-1][j])
      elif i-1 >= 0 and j-1 >= 0:
        distance_arr[i][j] = min(distance_arr[i][j],1 + distance_arr[i-1][j], 1 + distance_arr[i][j-1])

  for i in range(x_len-1, -1, -1):
    for j in range(y_len-1, -1, -1):
      if i+1 > x_len-1 and j+1 <=y_len-1:
        distance_arr[i][j] = min(distance_arr[i][j],1+distance_arr[i][j+1])
      elif i+1 <= x_len-1 and j+1 > y_len-1:
        distance_arr[i][j] = min(distance_arr[i][j],1+distance_arr[i+1][j])
      elif i+1 <= x_len-1 and j+1 <= y_len-1:
        distance_arr[i][j] = min(distance_arr[i][j],1+distance_arr[i+1][j], 1+distance_arr[i][j+1])

  return distance_arr
