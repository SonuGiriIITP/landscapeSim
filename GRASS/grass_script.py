import DEM_ascii_generator
import grass.script as g
import sys
import os

def main():
    resolution = 30
    DEM_generator(resolution)
    g.run_command('g.region', flags = 'ap', res = resolution)
    #Convert ASCII raster file to binary raster map layer
    pathname = os.path.dirname(sys.argv[0])        
    fullpath = os.path.abspath(pathname)
    g.run_command('r.in.ascii', overwrite = True, flags='f', input = fullpath + '/DEM.asc', output='FM2D_DEM_raw')
    g.run_command('r.out.png', input='FM2D_DEM_raw@user1', output = fullpath + '/'+'DEM')
    #Flow computation for massive grids (float version) 
    g.run_command('r.terraflow', overwrite = True, elevation = 'FM2D_DEM_raw@user1', filled = 'flooded_DEM',\
            direction = 'DEM_flow_direction' ,swatershed = 'DEM_sink_watershed' , accumulation = 'DEM_flow_accum' , tci = 'DEM_tci')
    g.run_command('r.out.png', input='flooded_DEM@user1', output = fullpath + '/'+'flooded')
    g.run_command('r.out.png', input='DEM_flow_direction@user1', output = fullpath + '/'+'direction')
    g.run_command('r.out.png', input='DEM_sink_watershed@user1', output = fullpath + '/'+'watershed')
    g.run_command('r.out.png', input='DEM_flow_accum@user1', output = fullpath + '/'+'accum')
    g.run_command('r.out.png', input='DEM_tci@user1', output = fullpath + '/'+'tci')
    #Generates raster maps of slope, aspect, curvatures and partial derivatives from a elevation raster map.
    #Aspect is calculated counterclockwise from east. 
    g.run_command('r.slope.aspect',overwrite=True,elevation='FM2D_DEM_raw@user1',slope='DEM_Slope',aspect='DEM_Aspect',dx='DEM_dx',dy='DEM_dy')
    g.run_command('r.out.png', input='DEM_Slope@user1', output = fullpath + '/'+'slope')
    g.run_command('r.out.png', input='DEM_Aspect@user1', output = fullpath + '/'+'aspect')
    g.run_command('r.out.png', input='DEM_dx@user1', output = fullpath + '/'+'dx')
    g.run_command('r.out.png', input='DEM_dy@user1', output = fullpath + '/'+'dy')
    #Generate images with textural features from a raster map. 
    g.run_command('r.texture', overwrite=True, flags = 'k',input = 'FM2D_DEM_raw@user1',prefix = 'texture',size = 5 )


def DEM_generator(resolution):
    #Parametes for genarating DEM 
    max_level = 10
    sigma = 1
    H = 0.85
    addition = True
    wrap = False
    gradient = 0
    seed = 31
    normalise=True
    bounds = [660,1840]
    gradient_values = [34,67,31,56]

    #Use FM2D algorithm to generate the DEM with specified parameters
    DEM_arr = DEM_ascii_generator.midPointFm2d(max_level,sigma,H,addition,wrap,gradient,seed,normalise,bounds,gradient_values)

    #write DEM to an ascii file
    fo = open("DEM.asc","w")
    north = 4928010
    south = north - resolution*DEM_arr.shape[0]
    east =  609000
    west =  east - resolution*DEM_arr.shape[0]
    rows = DEM_arr.shape[0]
    cols = DEM_arr.shape[1]
    null = -9999
    fo.write("north: "+str(north)+"\n")
    fo.write("south: "+str(south)+"\n")
    fo.write("east: "+str(east)+"\n")
    fo.write("west: "+str(west)+"\n")
    fo.write("rows: "+str(rows)+"\n")      
    fo.write("cols: "+str(cols)+"\n")
    fo.write("null: "+str(null)+"\n")
    for i in range(0,DEM_arr.shape[0]):
        fo.write("\n")
        for j in range(DEM_arr.shape[1]):
            fo.write(str(DEM_arr[i][j])+ " ")
    fo.close()

if __name__ == "__main__":
    main()
