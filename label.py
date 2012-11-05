from scipy.ndimage.measurements import (label, find_objects)
import pylab
import numpy

Agri_arr = numpy.load("Agri_arr.npy")
fields = numpy.load("Display_fields.npy")
s= [[1,1,1],
    [1,1,1],
    [1,1,1]]
mask = [Agri_arr == 1]
Agri_arr[mask] = 0
label_arr,no_of_patches =  label(Agri_arr,structure=s)
pylab.imsave("label_arr.png",label_arr,cmap="gray")
print "no of patches", no_of_patches
print numpy.max(Agri_arr)

