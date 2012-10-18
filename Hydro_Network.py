import numpy
import scipy
import math
import pylab
import MapGeneration_pure_python
import SpectralSynthesisFM2D
from scipy import ndimage
import city_block_dist_erosion

def pixel_exist(p,q,x_len,y_len):
    """
    Checks weather a pixel index is valid or not. 
    Args:
        p,q          : Index of the pixel to be checked (int, int)
        x_len, y_len : No of rows and columns in 2-D grid (int, int)
    Returns:
       True/False depending upon weather the 3x3 window centered at pixel(p,q) 
       lie within the boundary of the grid (boolean)
    """
    return ( (p >= 1) and (q >= 1) and ( p <= x_len - 2) and (q <= x_len-2) )


def DEM_creator(H, H_wt, seed, elev_range, max_level, DEMcreator_option):
    """
    Generates a DEM map with specified parameters described below.
    Args:
        H    : List containing H values         --> [H1, H2, H3, ...]( List of floats )
        H_wt : List containing H_wt             --> [H1_wt, H2_wt, H3_wt, ...] ( List of floats )
        seed : seed for random no generator     --> [seed1, seed2 seed3, ...] ( list of ints)
        elev_range: List containing elev bounds --> [elev_min, elev_max] (int, int) 
        max_level : size if grid 2^(max_level) + 1   ( int )
        DEMcreator_option: Specify the method used to generate DEM (fm2D or SS)       
    Result:
        Output : [ DEM_arr, [List of DEM input grids], [ List of DEM input grid names ], [list of parameters for DEM input grid]]
               where DEM_arr :Final Output DEM grid  
    """
    Input_DEMs = []
    name = []
    parameter = []
    if DEMcreator_option == 'fm2D':
        #Generate first DEM with gradient = 1 (i.e. TRUE) using FM2D method
        DEM_arr = MapGeneration_pure_python.midPointFm2d(max_level = max_level, sigma = 1, H = H[0], addition = True,\
                                          wrap = False, gradient = 1, seed = seed[0], normalise=True, bounds = elev_range)
        Input_DEMs.append(DEM_arr)
        file_name = "InputDEM_arr%d" % (1)
        name.append(file_name)
        parameter.append((H[0], H_wt[0]))
        DEM_arr = DEM_arr*H_wt[0]        
        for i in range(1,len(H)):    
            #Generate other DEM's with gradient = 0 (i.e. FLASE) and method specified by DEMcreator_option         
            temp_arr = MapGeneration_pure_python.midPointFm2d(max_level = max_level, sigma = 1, H = H[i], addition = True,\
                                                   wrap = False, gradient = 0, seed = seed[i], normalise = True, bounds = elev_range)
            Input_DEMs.append(temp_arr)
            file_name = "InputDEM_arr%d" % (i+1)
            name.append(file_name)
            parameter.append((H[i], H_wt[i]))
            DEM_arr = DEM_arr +  H_wt[i]*temp_arr     
    else:
        #Generate first DEM with gradient = 1 (i.e. TRUE) using FM2D method
        DEM_arr = MapGeneration_pure_python.midPointFm2d(max_level = max_level, sigma = 1, H = H[0], addition = True,\
                                          wrap = False, gradient = 1, seed = seed[0], normalise=True, bounds = elev_range)
        DEM_arr = DEM_arr[:-1,:-1] # omit last row and last col to maintain consistency with Spectral method output grid 
        Input_DEMs.append(DEM_arr)
        file_name = "InputDEM_arr%d" % (1)
        name.append(file_name)
        parameter.append((H[0], H_wt[0]))
        DEM_arr = H_wt[0]*DEM_arr

        for i in range(1,len(H)):
            #Generate other DEM's with gradient = 0 (i.e. FLASE) and method specified by DEMcreator_option i.e SS in this case    
            temp_arr = SpectralSynthesisFM2D.SpectralSynthesisFM2D(max_level = max_level, sigma = 1, H = H[i],\
                                                   seed = seed[i], normalise = True, bounds = elev_range)
            Input_DEMs.append(temp_arr)
            file_name = "InputDEM_arr%d" % (i+1)
            name.append(file_name)
            parameter.append((H[i], H_wt[i]))
            DEM_arr = DEM_arr +  H_wt[i]*temp_arr

    Output = [DEM_arr, Input_DEMs, name, parameter] 
    return Output


def Single_Cell_PitRemove(originalDEM, no_of_itr):
    """ 
    Removes single cell pits using 3x3 window
    Args:
        originalDEM: Digital Elevation Model (DEM) array ( 2-D array of floats )
        no_of_itr  : No of iteration of pit removal to be performed (int)
    Returns:
        OriginalDEM :Digital Elevation Model (DEM) with some pits removed using 3x3 window  ( 2-D array of floats )
    """
    #Get the dimensions of originalDEM
    (x_len,y_len) = originalDEM.shape
    for no_of_iteration in range(0, no_of_itr):
        #Run pit removal algorithm for no_of_itr times 
        for i in range(1,x_len-1):
            for j in range(1,y_len-1):
                #Remove pits using 3 x 3 window
                #Get the minimum of all values except central element
                A = min(originalDEM[i-1][j-1],originalDEM[i-1][j],originalDEM[i-1][j+1], originalDEM[i][j-1],\
                       originalDEM[i][j+1],originalDEM[i+1][j-1],originalDEM[i+1][j],originalDEM[i+1][j+1])
                if originalDEM[i][j] < A:
                    #If central pixel's value, in 3x3 window, is less then other pixels value
                    #in 3x3 window it's a pit.. Increase it's value to 'one + minimum of other 
                    # value i.e A'
                    originalDEM[i][j] = A + 1
    return originalDEM


def Flow_Dirn_3x3(DEM_arr, Flow_arr, pit_list):
    """
    Carry out flow direction calculation on DEM_arr using 3x3 window.
    Args:
        DEM_arr: Digital elevation array ( 2-D array of floats )
        Flow_arr: 2-D array of size DEM_arr to hold flow dirn ( 2-D array of ints)
        pit_list: An empty list required to hold pits ( Empty list [] )
    Result:
        pit_list: list containing single cell pits, will be used for catchment extraction ( list of tuples (int, int) )
        Flow_arr: 2-D array containing flow direction, will be used for catchment extraction (2-D array of ints )
    """
    pit_list = []
    (x_len, y_len) = DEM_arr.shape
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
    return (pit_list, Flow_arr)


def Erosion(River_arr, DEM_arr, river_drop):
    """ 
    Carry out erosion based on city-block Distance from River 
    Args:
        River_arr : 2-D array containing flow accumulation values ( 2-D array of ints )
        DEM_arr   : DEM array (an array of floats)
        river_drop: describes the maximum extent of erosion to be performed
                    A proportion of this value is subtracted from pixel depending on distance from the river (float)
    Result:
        DEM_arr: Eroded digital elevation model(DEM) array ( 2-D array of floats )
        Distance_arr: 2-D array having nearest city-block distance of a pixel form River (2-D array of int)
    """
    (x_len, y_len) = DEM_arr.shape
    Distance_arr = city_block_dist_erosion.CityBlock(River_arr)
    # Create a mask for differnet distances used for DEM erosion
    mask4 = [Distance_arr <= 15]
    mask5 = [Distance_arr > 3]
    mask3 = [Distance_arr == 3]
    mask2 = [Distance_arr == 2]
    mask1 = [Distance_arr == 1]
    mask0 = [Distance_arr == 0]
    max_flow_accum = numpy.max(River_arr)
    for i in range(0,x_len):
        for j in range(0,y_len):
        #Erode the landscape based on it's distance from river
            if mask0[0][i][j] == True:
                DEM_arr[i][j] = DEM_arr[i][j] - river_drop
            elif mask1[0][i][j] == True:
                DEM_arr[i][j] = DEM_arr[i][j] - (river_drop * 0.8)
            elif mask2[0][i][j] == True:
                DEM_arr[i][j] = DEM_arr[i][j] - (river_drop * 0.6)
            elif mask3[0][i][j] == True:
                DEM_arr[i][j] = DEM_arr[i][j] - (river_drop * 0.4)
            elif mask4[0][i][j] == True and mask5[0][i][j] == True:
                DEM_arr[i][j] = DEM_arr[i][j] - (river_drop * 0.2)
            else:     
                DEM_arr[i][j] = DEM_arr[i][j] - (river_drop * 0.1)
    return (DEM_arr, Distance_arr)


def CatchmentExtraction(pit_list, DEM_arr, Flow_arr, max_posn):
    """
    Carry out catchment extraction and depression fillings
    Args:
        pit_list: list containing pit positions identified using 3x3 window (empty list [])  
        DEM_arr : DEM array neede for depression filling ( 2-D array of floats )
        Flow_arr: Flow array containing flow directions ( 2-D array of int )
        max_posn: position of pixel having highest elevation in DEM_arr ( tuple (int, int) )
    Result:
        (DEM_arr, Found_arr, Catchment_boundary_arr)
        DEM_arr: DEM_arr with some depressions filled (2-D array of flaot)
        Found_arr: Actual catchment array (2-D arrays of ints)
        Catchment_boundary_arr: Catchment boundary(inner) arrays(2-D array of ints)
    """
    # label will be used to assign labels to differnet catchments
    label = 0 
    (x_len,y_len) = DEM_arr.shape
    # Found_arr will hold the catchment with different labels
    Found_arr = numpy.zeros( (x_len,y_len),dtype = "uint8")
    # Catchment_boundary_arr will hold the Catchment boundaries
    Catchment_boundary_arr = numpy.zeros((x_len,y_len) , dtype = "uint8" )
    while len(pit_list) >= 1:
    #___________For each and every pit in the DEM do ___________
        stack = []
        pit = pit_list.pop(0)
        stack.append(pit)
        #Increase the label being assigned to the catchment_______
        label = label + 1  
        #_______Identify catchment for each and every pit_________
        catchment_pixels = []
        catchment_pixels.append(( DEM_arr[pit[0],pit[1]], pit[0], pit[1] ))
        while len(stack) > 0:
            (p,q) = stack.pop(0)
            Found_arr[p][q] = label
            # Pop an element from stack check if its adjacent pixels exist and contribute 
            # its flow to the central pixel(pixel popped) then append it into list, continue
            # this untill stack gets empty
            if pixel_exist(p-1,q-1,x_len,y_len):
                if Flow_arr[p-1][q-1] == 7:
                    catchment_pixels.append((DEM_arr[p-1,q-1],p-1,q-1))
                    stack.append((p-1,q-1))
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
        pour_point = max_posn
        flag = 0
        for i in range(0,len(catchment_pixels)):
            (p,q) = ( catchment_pixels[i][1],catchment_pixels[i][2] )
            label = Found_arr[p][q]
            # Catchment Outlet(pour point) will be the minimum catchment boundary pixel
            if (Found_arr[p-1][q-1] != label or Found_arr[p-1][q] != label or Found_arr[p-1][q+1] != label or 
                Found_arr[p][q-1] != label or Found_arr[p][q+1] != label or Found_arr[p+1][q-1] != label or
                Found_arr[p+1][q] != label or Found_arr[p+1][q+1] != label):# if pixel lie on boundary of catchment
                Catchment_boundary_arr[p][q] = 255 #Inner boundary gets marked
                if DEM_arr[ pour_point[0] ][ pour_point[1] ] > DEM_arr[p][q]:#if height of boundary is less then update pour point
                    pour_point = (p,q)
                    flag = 1
        if flag == 1:
            #Pour_point_list.append((DEM_arr[pour_point],pour_point[0],pour_point[1]))
            #Pour_point_arr[pour_point] = 255
            for i in range(0,len(catchment_pixels)):
                if catchment_pixels[i][0] < DEM_arr[pour_point]:
                #Fill the depression in the catchment
                    DEM_arr[catchment_pixels[i][1],catchment_pixels[i][2]] = DEM_arr[pour_point]

    return (DEM_arr, Found_arr, Catchment_boundary_arr)


def Get_Flow_Dirn_using_9x9_window(DEM, Flow_dirn_arr , pit_list):
    """
    Given a DEM the function returns the flow direction matrix, it also erodes the DEM
    as per requirement to direct flow  
    Args:
        DEM:  Digital Elevation Model (2-D array of floats)
        Flow_dirn_arr: Empty flow direction array having all entries zero (2-D array of tuples (0,0) )
        pit_list :Empty list used to hold pits ( Empty List )
    Result:
        pit_list: pits found using 9x9 window ( List of tuple (int, int) )
        Flow_dirn_arr: 2-D array containing Flow Directions  (2-D array of tuples (int, int) )
        DEM: Modified DEM after little erosion during flow direction assignment (2-D array of floats)
    """
    (x_len,y_len) = DEM.shape
    pit_list = []
  
    #Get the flow direction using 9x9 window
    for i in range(4,x_len-4):
        for j in range(4,y_len-4):#loop index start from 4 and ends at len-4 to handle boundary cases
            if Flow_dirn_arr[i][j][0] == 0 and Flow_dirn_arr[i][j][1] == 0:
                (x,y) = ndimage.minimum_position( DEM[i - 4:i + 5, j - 4:j + 5] )
                # (x,y) is the position of minimum element in 9x9 window
                (min_x,min_y) =  (x - 4, y - 4)
                # (min_x,min_y) is the position of minimum element in 9x9 window with origin
                # shifted to the central pixel

                # 9x9 window can be divided into 4 quadrants ,7 lines of code below takes care of 3 
                # other quadrants in 9x9 window, since we are writing a general code for first quadrant,
                # where q and p are non-negative integers 
                sign_x = 1 # indicative of +ve x value
                sign_y = 1 # indicative of +ve y value
                if min_x < 0:
                    sign_x = -1
                if min_y < 0:
                    sign_y = -1
                (p, q) = (abs(min_x),abs(min_y)) 

                Elev_diff = (DEM[i][j] - DEM[p*sign_x + i][q*sign_y + j])/max(p,q)
                # difference in elevation of the central pixel and the pixel with minimum elevation
                # in 9x9 window, required for the purpose of erosion 
     
                #Different cases in the 9x9 window has been handled in various if-else statements
                if p == 0:
                    if q == 1:
                        Flow_dirn_arr[i][j] = (i + 0*sign_x,j + 1*sign_y)
                    elif q == 2:
                        Flow_dirn_arr[i][j] = (i + 0*sign_x,j + 1*sign_y)
                        Flow_dirn_arr[i + 0*sign_x][j + 1*sign_y] = (0*sign_x + i ,2*sign_y + j)
                        DEM[i + 0*sign_x][j + 1*sign_y] = DEM[0*sign_x + i][2*sign_y + j] + Elev_diff            
                    elif q == 3:
                        Flow_dirn_arr[i][j] = (i + 0*sign_x,j + 1*sign_y)
                        Flow_dirn_arr[i + 0*sign_x][j + 1*sign_y] = (0*sign_x + i ,2*sign_y + j)
                        Flow_dirn_arr[i + 0*sign_x][j + 2*sign_y] = (0*sign_x + i, 3*sign_y + j)
                        DEM[i + 0*sign_x][j + 1*sign_y] = DEM[0*sign_x + i][3*sign_y + j] + 2*Elev_diff
                        DEM[i + 0*sign_x][j + 2*sign_y] = DEM[0*sign_x + i][3*sign_y + j] + Elev_diff
                    elif q == 4:
                        Flow_dirn_arr[i][j] = (i + 0*sign_x,j + 1*sign_y)
                        Flow_dirn_arr[i + 0*sign_x][j + 1*sign_y] = (0*sign_x + i ,2*sign_y + j)
                        Flow_dirn_arr[i + 0*sign_x][j + 2*sign_y] = (0*sign_x + i, 3*sign_y + j)
                        Flow_dirn_arr[i + 0*sign_x][j + 3*sign_y] = (0*sign_x + i, 4*sign_y + j)
                        DEM[i + 0*sign_x][j + 1*sign_y] = DEM[0*sign_x + i][4*sign_y + j] + 3*Elev_diff
                        DEM[i + 0*sign_x][j + 2*sign_y] = DEM[0*sign_x + i][4*sign_y + j] + 2*Elev_diff
                        DEM[i + 0*sign_x][j + 3*sign_y] = DEM[0*sign_x + i][4*sign_y + j] + Elev_diff

                if p == 1:
                    if q == 0:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 0*sign_y )
                    elif q == 1:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                    elif q == 2:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 1*sign_x ,j + 2*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 1*sign_x][j + 2*sign_y] + Elev_diff
                    elif q == 3:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 1*sign_x ,j + 2*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 2*sign_y] = (i + 1*sign_x ,j + 3*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 1*sign_x ][j + 3*sign_y ] + 2*Elev_diff
                        DEM[i + 1*sign_x][j + 2*sign_y] = DEM[i + 1*sign_x ][j + 3*sign_y ] + Elev_diff
                    elif q == 4:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 1*sign_x ,j + 2*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 2*sign_y] = (i + 1*sign_x ,j + 3*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 3*sign_y] = (i + 1*sign_x ,j + 4*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 1*sign_x ][j + 4*sign_y ] + 3*Elev_diff
                        DEM[i + 1*sign_x][j + 2*sign_y] = DEM[i + 1*sign_x ][j + 4*sign_y ] + 2*Elev_diff
                        DEM[i + 1*sign_x][j + 3*sign_y] = DEM[i + 1*sign_x ][j + 4*sign_y ] + Elev_diff     

                if p == 2:
                    if q == 0:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 0*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 0*sign_y] = (i + 2*sign_x ,j + 0*sign_y )
                        DEM[i + 1*sign_x][j + 0*sign_y] = DEM[i + 2*sign_x ][j + 0*sign_y] + Elev_diff  
                    elif q == 1:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 1*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 2*sign_x ][j + 1*sign_y] + Elev_diff
                    elif q == 2:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 2*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 2*sign_x ][j + 2*sign_y ] + Elev_diff
                    elif q == 3:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 2*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 2*sign_y] = (i + 2*sign_x ,j + 3*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 2*sign_x ][j + 3*sign_y ] + 2*Elev_diff
                        DEM[i + 2*sign_x][j + 2*sign_y] = DEM[i + 2*sign_x ][j + 3*sign_y ] + Elev_diff
                    elif q == 4:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 2*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 2*sign_y] = (i + 2*sign_x ,j + 3*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 3*sign_y] = (i + 2*sign_x ,j + 4*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 2*sign_x ][j + 4*sign_y ] + 3*Elev_diff
                        DEM[i + 2*sign_x][j + 2*sign_y] = DEM[i + 2*sign_x ][j + 4*sign_y ] + 2*Elev_diff
                        DEM[i + 2*sign_x][j + 3*sign_y] = DEM[i + 2*sign_x ][j + 4*sign_y ] + Elev_diff

                if p == 3:
                    if q == 0:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 0*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 0*sign_y] = (i + 2*sign_x ,j + 0*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 0*sign_y] = (i + 3*sign_x ,j + 0*sign_y )
                        DEM[i + 1*sign_x][j + 0*sign_y] = DEM[i + 3*sign_x][j + 0*sign_y] + 2*Elev_diff
                        DEM[i + 2*sign_x][j + 0*sign_y] = DEM[i + 3*sign_x][j + 0*sign_y] + Elev_diff 
                    elif q == 1:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 1*sign_y] = (i + 3*sign_x ,j + 1*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 3*sign_x][j + 1*sign_y] + 2*Elev_diff
                        DEM[i + 2*sign_x][j + 1*sign_y] = DEM[i + 3*sign_x][j + 1*sign_y] + Elev_diff
                    elif q == 2:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 2*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 2*sign_y] = (i + 3*sign_x ,j + 2*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 3*sign_x][j + 2*sign_y] + 2*Elev_diff
                        DEM[i + 2*sign_x][j + 2*sign_y] = DEM[i + 3*sign_x][j + 2*sign_y] + Elev_diff
                    elif q == 3:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 2*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 2*sign_y] = (i + 3*sign_x ,j + 3*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 3*sign_x ][j + 3*sign_y ] + 2*Elev_diff
                        DEM[i + 2*sign_x][j + 2*sign_y] = DEM[i + 3*sign_x ][j + 3*sign_y ] + Elev_diff
                    elif q == 4:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 2*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 2*sign_y] = (i + 3*sign_x ,j + 3*sign_y )
                        Flow_dirn_arr[i + 3*sign_x][j + 3*sign_y] = (i + 3*sign_x ,j + 4*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 3*sign_x ][j + 4*sign_y ] + 3*Elev_diff
                        DEM[i + 2*sign_x][j + 2*sign_y] = DEM[i + 3*sign_x ][j + 4*sign_y ] + 2*Elev_diff
                        DEM[i + 3*sign_x][j + 3*sign_y] = DEM[i + 3*sign_x ][j + 4*sign_y ] + Elev_diff
     
                if p == 4:
                    if q == 0:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 0*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 0*sign_y] = (i + 2*sign_x ,j + 0*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 0*sign_y] = (i + 3*sign_x ,j + 0*sign_y ) 
                        Flow_dirn_arr[i + 3*sign_x][j + 0*sign_y] = (i + 4*sign_x ,j + 0*sign_y )
                        DEM[i + 1*sign_x][j + 0*sign_y] = DEM[i + 4*sign_x ][j + 0*sign_y ] + 3*Elev_diff
                        DEM[i + 2*sign_x][j + 0*sign_y] = DEM[i + 4*sign_x ][j + 0*sign_y ] + 2*Elev_diff
                        DEM[i + 3*sign_x][j + 0*sign_y] = DEM[i + 4*sign_x ][j + 0*sign_y ] + Elev_diff
                    elif q == 1:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 1*sign_y] = (i + 3*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 3*sign_x][j + 1*sign_y] = (i + 4*sign_x ,j + 1*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 4*sign_x ][j + 1*sign_y ] + 3*Elev_diff
                        DEM[i + 2*sign_x][j + 1*sign_y] = DEM[i + 4*sign_x ][j + 1*sign_y ] + 2*Elev_diff
                        DEM[i + 3*sign_x][j + 1*sign_y] = DEM[i + 4*sign_x ][j + 1*sign_y ] + Elev_diff
                    elif q == 2:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 2*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 2*sign_y] = (i + 3*sign_x ,j + 2*sign_y )
                        Flow_dirn_arr[i + 3*sign_x][j + 2*sign_y] = (i + 4*sign_x ,j + 2*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 4*sign_x ][j + 2*sign_y ] + 3*Elev_diff
                        DEM[i + 2*sign_x][j + 2*sign_y] = DEM[i + 4*sign_x ][j + 2*sign_y ] + 2*Elev_diff
                        DEM[i + 3*sign_x][j + 2*sign_y] = DEM[i + 4*sign_x ][j + 2*sign_y ] + Elev_diff
                    elif q == 3:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 2*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 2*sign_y] = (i + 3*sign_x ,j + 3*sign_y )
                        Flow_dirn_arr[i + 3*sign_x][j + 3*sign_y] = (i + 4*sign_x ,j + 3*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 4*sign_x ][j + 3*sign_y ] + 3*Elev_diff
                        DEM[i + 2*sign_x][j + 2*sign_y] = DEM[i + 4*sign_x ][j + 3*sign_y ] + 2*Elev_diff
                        DEM[i + 3*sign_x][j + 3*sign_y] = DEM[i + 4*sign_x ][j + 3*sign_y ] + Elev_diff
                    elif q == 4:
                        Flow_dirn_arr[i][j] = (i + 1*sign_x ,j + 1*sign_y )
                        Flow_dirn_arr[i + 1*sign_x][j + 1*sign_y] = (i + 2*sign_x ,j + 2*sign_y )
                        Flow_dirn_arr[i + 2*sign_x][j + 2*sign_y] = (i + 3*sign_x ,j + 3*sign_y )
                        Flow_dirn_arr[i + 3*sign_x][j + 3*sign_y] = (i + 4*sign_x ,j + 4*sign_y )
                        DEM[i + 1*sign_x][j + 1*sign_y] = DEM[i + 4*sign_x ][j + 4*sign_y ] + 3*Elev_diff
                        DEM[i + 2*sign_x][j + 2*sign_y] = DEM[i + 4*sign_x ][j + 4*sign_y ] + 2*Elev_diff
                        DEM[i + 3*sign_x][j + 3*sign_y] = DEM[i + 4*sign_x ][j + 4*sign_y ] + Elev_diff
                if p == 0 and q == 0:
                    pit_list.append((i,j)) 
    return (pit_list, Flow_dirn_arr,DEM)


def Flow_accumulation(Flow_dirn_arr, Flow_Accum_arr, DEM):
    """
    Performs flow accumulation over DEM using Flow_dirn_arr generated using 9x9 window
    Args:
        Flow_dirn_arr: Array containing flow directions, each cell contains the 
                       index ((i,j) values) of the pixel where flow is to be directed (2-D array of tuple of the form(int,int))
        Flow_Accum_arr: Array to hold flow accumulation values (2-D array of ints)
        DEM: Eroded Digital Elevation Model (2-D array of floats)
    Result:
        Flow_Accum_arr : Matrix containing flow accumulation values (2-D array of ints)
    """
    (x_len, y_len, z) = Flow_dirn_arr.shape
    # Create and initialze a list containing tuples (DEM[i][j],i,j)
    A = []
    for i in range(0,x_len):
        for j in range(0,y_len): 
            A.append( (DEM[i][j],i,j) )
  
    # Do flow accumulation
    A.sort(reverse = 1 )
    for i in range( 0, len(A) ):
        #for each pixel in decreasing order of height...
        if Flow_dirn_arr[ A[i][1] , A[i][2], 0 ] != 0 and Flow_dirn_arr[ A[i][1] , A[i][2], 1 ] != 0:
            (x,y) = (Flow_dirn_arr[ A[i][1] , A[i][2]])
            Flow_Accum_arr[x][y] = Flow_Accum_arr[x][y] + Flow_Accum_arr[ A[i][1] ][ A[i][2] ]
    return Flow_Accum_arr
