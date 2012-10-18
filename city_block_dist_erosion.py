import numpy
from scipy import *

def CityBlock(image_arr):
#converting the gray-scale image into binary image
  (x_len,y_len)=image_arr.shape
  L = []
  for i in range(0,x_len):
    for j in range(0,y_len):
      L.append((image_arr[i][j],i,j))

  L.sort()
  binary_arr = numpy.zeros((x_len, y_len),dtype="uint8") 
  for i in range(len(L)-1 , 95*len(L)/100,-1):
        binary_arr[ L[i][1] ][ L[i][2] ] = 255

#Initializing distance array  
  distance_arr=zeros((x_len,y_len), dtype= "int")
  for i in range(0, x_len):
    for j in range(0, y_len):
      if binary_arr[i][j]==0:
        distance_arr[i][j] = x_len + y_len+1
      else:
        distance_arr[i][j] = 0

#Applying distance transform
  for i in range(0, x_len):
    for j in range(0, y_len):
      if i-1 < 0 and j-1 >=0:
        distance_arr[i][j] = min(distance_arr[i][j],1+distance_arr[i][j-1]);
      elif i-1 >= 0 and j-1 < 0:
        distance_arr[i][j] = min(distance_arr[i][j],1+distance_arr[i-1][j]);
      elif i-1 >= 0 and j-1 >= 0:
        distance_arr[i][j] = min(distance_arr[i][j],1+distance_arr[i-1][j], 1+distance_arr[i][j-1]);
      else:
        distance_arr[i][j]=distance_arr[i][j];

  for i in range(x_len-1, -1, -1):
    for j in range(y_len-1, -1, -1):
      if i+1 > x_len-1 and j+1 <=y_len-1:
        distance_arr[i][j] = min(distance_arr[i][j],1+distance_arr[i][j+1]);
      elif i+1 <= x_len-1 and j+1 > y_len-1:
        distance_arr[i][j] = min(distance_arr[i][j],1+distance_arr[i+1][j]);
      elif i+1 <= x_len-1 and j+1 <= y_len-1:
        distance_arr[i][j] = min(distance_arr[i][j],1+distance_arr[i+1][j], 1+distance_arr[i][j+1]);
      else:
        distance_arr[i][j]=distance_arr[i][j];
  return distance_arr
