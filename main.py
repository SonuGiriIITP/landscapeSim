#!/usr/bin/env python
import pylab
import numpy
from scipy import ndimage 
import time
import yaml
import Hydro_Network
import DecisionTree
import VegetationClassify 
import Geometry
import surface_plot

def main():
    """
    It imports all necessary python modules required for landscape simulation.
    It then performs the following:
    1. Gets all the parameters required for simulation from parameter.yaml file. 
    2. calls DEM_creator() --> for generating DEM grid
    3. Iteratively operate on DEM grid and do the following
       3.1 Remove single cell pits by calling Single_Cell_PitRemove()
       3.2 Get flow dirn using 9x9 window by calling Get_Flow_Dirn_using_9x9_window()  
       3.3 Get flow dirn using 3x3 window by calling Flow_Dirn_3x3() for catchment extraction
       3.4 Extract the catchment and do depression filling using CatchmentExtraction()
       3.5 Again get flow direction using Get_Flow_Dirn_using_9x9_window() after depression filling
       3.6 Perfrom flow accumulation by calling Flow_accumulation() 
       3.7 Do the erosion by calling Erosion()
    4. Generate a Decision tree for land_cover allocation by calling DecisionTree()
    5. Assign the vegetation class to DEM by calling VegetationClassify()
    6. Generate some agricultural field by calling GeometricFeature()
    """
    time1 = time.time()
    # Get the parameters for parameter.yaml file
    yaml_file = open('Parameters/parameters.yaml', 'r')
    stream = yaml.load(yaml_file)
    counter = stream['counter']
    H = stream['H']
    H_wt =  stream['H_wt']
    seed = stream['seed']
    elev_range = stream['elev_range']
    river_drop = stream['river_drop']
    max_level = stream['max_level']
    DEMcreator_option = stream['DEMcreator_option']
    Three_DplotDEM = stream["Three_DplotDEM"]
    output_dir = stream['output_dir']
    response = stream['response']
    min_area = stream['min_area']
    max_area = stream['max_area']
    aspect_ratio = stream['aspect_ratio']
    agri_area_limit = stream['agri_area_limit']
    print ("Running simulation with follwing parameters")
    print ("Counter %d" % counter)
    print ("H %s" % H)
    print ("H_wt %s" % H_wt)
    print ("seed %s" % seed)
    print ("elev_range %s" % elev_range)
    print ("river_drop %s" % river_drop)
    print ("max_level %s" % max_level)
    print ("DEMcreator_option %s" % DEMcreator_option)
    print ("output_dir %s" % output_dir)
    print ("response %s" % response)
    print ("min_area %d" % min_area)
    print ("max_area %d" % max_area)
    print ("aspect_ratio %s" % aspect_ratio)
    print ("agri_area_limit %s" % agri_area_limit)

    #Generate DEM using FM2D/SS algorithm by calling DEM_creator(args...) function")
    DEM_Result = Hydro_Network.DEM_creator(H, H_wt, seed, elev_range, max_level, DEMcreator_option)
    #Write result to Output file
    file_name = "%s/Original_DEM" % (output_dir)
    pylab.imsave(file_name, DEM_Result[0])
    for i in range(0,len(DEM_Result[1])):
        file_name = "%s/%s" % (output_dir,DEM_Result[2][i])#,DEM_Result[3][i][0],DEM_Result[3][i][1])
        pylab.imsave(file_name, DEM_Result[1][i])

    DEM = DEM_Result[0]
    for iteration in range(0,counter):
        #Remove sink using 3x3 window by calling Single_Cell_PitRemove(originalDEM, no_of_itr)
        DEM = Hydro_Network.Single_Cell_PitRemove(DEM, no_of_itr = 6)
        (x_len,y_len) = DEM.shape
        max_posn = ndimage.maximum_position(DEM)
        Flow_dirn_arr = numpy.zeros((x_len,y_len,2), dtype="int" )
        #Flow_arr will be used for the purpose of catchment extraction
        Flow_arr = numpy.zeros((x_len, y_len), dtype = "uint8")
        River_arr = numpy.ones((x_len, y_len), dtype = "int")
        pit_list = [] #Not required now
        ( pit_list, Flow_dirn_arr, DEM ) = Hydro_Network.Get_Flow_Dirn_using_9x9_window(DEM, Flow_dirn_arr, pit_list)
        # call Flow_Dirn_3x3(DEM, Flow_arr , pit_list) for the purpose of catchment extraction
        pit_list = [] #Required for catchment extraction
        ( pit_list, Flow_arr ) = Hydro_Network.Flow_Dirn_3x3(DEM, Flow_arr , pit_list) 
        #Catchment extraction, calling CatchmentExtraction(pit_list, DEM_arr, max_posn)
        (DEM, Found_arr, Catchment_boundary_arr) = Hydro_Network.CatchmentExtraction(pit_list, DEM, Flow_arr, max_posn)
        #Write result to Output file
        file_name = "%s/Catchment%s" % (output_dir, iteration+1)
        pylab.imsave(file_name, Found_arr)
        file_name = "%s/Catchment_Boundary%s" % (output_dir, iteration+1)
        pylab.imsave(file_name, Catchment_boundary_arr)        
        #Assignnig flow dirnection again after catchment extraction and Depression filling
        ( pit_list, Flow_dirn_arr, DEM ) = Hydro_Network.Get_Flow_Dirn_using_9x9_window(DEM, Flow_dirn_arr , pit_list)
        #Calculate flow accumulation by Calling Flow_accumulation(Flow_dirn_arr ,River_arr , DEM)
        River_arr = Hydro_Network.Flow_accumulation(Flow_dirn_arr ,River_arr, DEM)
        #Write result to Output file
        file_name = "%s/River%s" % (output_dir,iteration+1)
        pylab.imsave(file_name, River_arr)
        #"Eroding the DEM based on Distance form River ...Calling Erosion(River_arr,DEM_arr,river_drop)
        (DEM, Distance_arr) = Hydro_Network.Erosion(River_arr, DEM, river_drop)  
        #Write result to Output file
        file_name = "%s/ErodedDEM%s" % (output_dir, iteration+1)
        pylab.imsave(file_name, DEM)
        file_name = "%s/RiverDistance%s" % (output_dir, iteration+1)
        pylab.imsave(file_name, Distance_arr)
    
    if Three_DplotDEM == 'y' or Three_DplotDEM == 'Y':
        surface_plot.plot(DEM)
    time2 = time.time()
    print ("Time taken in Erosion modeling", time2 - time1,"seconds")

    if (response == 'y') or (response == 'Y'):
        elev_filename = stream["training_data_elev"]
        landcover_filename = stream["training_data_landcover"]
        river_filename = stream["training_data_river"]

        DecisionTree.DecisionTree(output_dir, elev_filename, landcover_filename, river_filename)
        time3 = time.time()
        print "Time taken to generate decision tree is " , (time3 - time2) ,"seconds"

    time3 = time.time()
    Veg_arr = VegetationClassify.VegetationClassify(DEM, River_arr)
    file_name = "%s/Landcover" % (output_dir)
    pylab.imsave(file_name, Veg_arr)
    time4 = time.time()
    print "Time taken to assign landcover is " , (time4 - time3),"seconds"
    fields =Geometry.GeometricFeature(Veg_arr, min_area = min_area, max_area = max_area, aspect_ratio = aspect_ratio ,agri_area_limit = agri_area_limit)
    file_name = "%s/Fields" % (output_dir)
    pylab.imsave(file_name, fields)
    time5 = time.time()
    print "Time taken to generate Geometric Features is " ,(time5 - time4) ,"seconds"    
    yaml_file.close()

if __name__ == "__main__":
    main()
