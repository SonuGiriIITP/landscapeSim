#!/usr/bin/python
#Import Packages (Note the following packages need to be present)
import scipy,numpy

def HistEqualization(image_arr):
  # Histogram array will contain the image_arr histogram
  histogram_arr = np.zeros(101)
  (x_len,y_len) = image_arr.shape
 
  for i in range(0, x_len):
    for j in range(0, y_len):
      # Calculating histogram from image_arr
      histogram_arr[int(image_arr[i][j])]=histogram_arr[int(image_arr[i][j])]+1
  
  # Compute cummulative histogram array
  cumm_histogram_arr = np.zeros(101)
  cumm_histogram_arr[0] = histogram_arr[0]
  for i in range(1, 101):
     cumm_histogram_arr[i] = cumm_histogram_arr[i-1]+ histogram_arr[i]

  # Computing the Look Up Table(LUT) from commulative histogram array
  LUT = np.zeros(101)
  factor = 101.0/(x_len*y_len)
  for i in range(101):
    LUT[i]= round(factor * cumm_histogram_arr[i])
  
  # Geting the new equalised array with the help of LUT  
  equalized_arr = numpy.zeros((x_len, y_len))
  for i in range(0, x_len):
    for j in range(0, y_len):
      equalized_arr[i][j] = LUT[int(image_arr[i][j])];

  return equalized_arr
