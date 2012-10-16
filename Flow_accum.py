import numpy
def Generate_River( Flow_arr,River_arr,DEM_arr):
  """Flow Accumulation step"""
  L = []
  for i in range(1,DEM_arr.shape[0]-1):
    for j in range(1,DEM_arr.shape[1]-1):  
      L.append((DEM_arr[i][j],i,j))
  L.sort()
 
  for p in range( len(L)-1 , -1 , -1 ):
    (value,i,j) = L[p]
    if Flow_arr[i][j] == 0:
      River_arr[i][j+1] = River_arr[i][j+1]  + River_arr[i][j]
    elif Flow_arr[i][j] == 1:
      River_arr[i-1][j+1] = River_arr[i-1][j+1] + River_arr[i][j] 
    elif Flow_arr[i][j] == 2:
      River_arr[i-1][j] = River_arr[i-1][j] + River_arr[i][j]
    elif Flow_arr[i][j] == 3:
      River_arr[i-1][j-1] = River_arr[i-1][j-1] + River_arr[i][j]
    elif Flow_arr[i][j] == 4:
      River_arr[i][j-1] = River_arr[i][j-1] + River_arr[i][j]
    elif Flow_arr[i][j] == 5:
      River_arr[i+1][j-1] = River_arr[i+1][j-1] + River_arr[i][j]          
    elif Flow_arr[i][j] == 6:
      River_arr[i+1][j] = River_arr[i+1][j] + River_arr[i][j]
    elif Flow_arr[i][j] == 7:
      River_arr[i+1][j+1] = River_arr[i+1][j+1] + River_arr[i][j]
    elif Flow_arr[i][j] == 8:
      elev = DEM_arr[i][j]
      River_arr[i][j] = numpy.max(River_arr[i-1:i+2,j-1:j+2])
      if River_arr[i-1][j-1] == elev:
        River_arr[i-1][j-1] = River_arr[i][j]
      elif River_arr[i-1][j] == elev:
        River_arr[i-1][j] = River_arr[i][j]  
      elif River_arr[i-1][j+1] == elev:
        River_arr[i-1][j+1] = River_arr[i][j]
      elif River_arr[i][j-1] == elev:
        River_arr[i][j-1] = River_arr[i][j]
      elif River_arr[i][j+1] == elev:
        River_arr[i][j+1] = River_arr[i][j]
      elif River_arr[i+1][j-1] == elev:
        River_arr[i+1][j-1] = River_arr[i][j]
      elif River_arr[i+1][j] == elev:
        River_arr[i+1][j] = River_arr[i][j]
      if River_arr[i+1][j+1] == elev:
        River_arr[i+1][j+1] = River_arr[i][j]
  return River_arr
