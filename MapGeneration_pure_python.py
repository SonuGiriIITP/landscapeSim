import numpy as np
import random

def fN(delta, x, sigma):
    '''Takes an n length numpy np.array, x, and 'blurs' the mean
       using a Gaussian distrubition centered about 0. The
       Gaussian distribution is scaled by delta and has standard
       deviation of sigma.'''
    return_num = x.mean() + delta * random.gauss(0.0, sigma)
    return return_num


def midPointFm2d(max_level, sigma, H, addition, wrap, gradient,
                 seed=0, normalise = True, bounds = [0,1]):
    """
    ________________________________________________________________________
    Args:
        max_level : Maximum number of recursions( N = 2^max_level)
        sigma     : Initial standard deviation
        H         : Roughness constant varies form 0.0 to 1.0
        addition  : boolean parameter (turns random additions on/off)
        wrap      : wraps the Image
        gradient  : if 1 then corners are deterministically set else randomally
        seed      : seed value for random number generator
        normalise : normalizes the data using bound
        bounds    : used for normalization of the grid data
    
    Result:     
        Output is given in the form of an array(grid) which holds surface
        elevation data for a square region.  
    _________________________________________________________________________
    """	
    
    random.seed(seed) #seed the random number generator 
    north      = 60
    north_west = 20
    west       = 5
    south_west = 25
    south      = 45
    south_east = 55
    east       = 25
    north_east = 60
    center     = 52

    N = 2**max_level 
    grid = np.zeros([N+1,N+1])#Generate a 2-D real array of size (N+1)^2 initialized to zero 

    delta = sigma # delta is a real variable holding standard deviations

    if gradient == 1: #set the initial random corners
        grid[0,0] = north_west
        grid[0,N/2] = north 
        grid[0,N] = north_east
        grid[N/2,N] = east
        grid[N,N] = south_east
        grid[N,N/2] = south
        grid[N,0] = south_west
        grid[N/2,0] = west
        grid[N/2,N/2] = center
        
    else:
        grid[0,0] = delta*random.gauss(0.0,sigma)
        grid[0,N/2] = delta*random.gauss(0.0,sigma) 
        grid[0,N] = delta*random.gauss(0.0,sigma)
        grid[N/2,N] = delta*random.gauss(0.0,sigma)
        grid[N,N] = delta*random.gauss(0.0,sigma)
        grid[N,N/2] = delta*random.gauss(0.0,sigma)
        grid[N,0] = delta*random.gauss(0.0,sigma)
        grid[N/2,0] = delta*random.gauss(0.0,sigma)
        grid[N/2,N/2] = delta*random.gauss(0.0,sigma)


    D = N/2 #D is maximum no of recursions (earlier it was N)
    dd = N/4 #modified code (earlier it was N/2)

    def locFN(x):
        return fN(delta, x, sigma)

    for stage in np.arange(2,max_level+1):#code modified here (earlier (1,max_level+1))
        
        delta = delta / 2**(H/2.0) #going from grid type I to grid Type II

        vec = range(dd, N-dd+1, D)
        #interpolate and offset points
        grid[dd:N-dd+1:D,dd:N-dd+1:D] = [[fN(delta, np.array([grid[x+dd,y+dd],\
			grid[x+dd,y-dd],grid[x-dd,y+dd], grid[x-dd,y-dd]]), sigma)\
            for y in vec] for x in vec]

        if addition:
            vec = range(0,N+1,D)
            grid[0:N+1:D,0:N+1:D] += [[delta*random.gauss(0.0,sigma)\
				for y in vec]for x in vec]

        delta = delta / 2**(H/2.0) #going from grid type II to grid type I

        for x in range(dd,N-dd+1,D):#interpolate and offset boundary grid points
            grid[x,0] = fN(delta, np.array([grid[x+dd,0],grid[x-dd,0],\
				grid[x,dd]]),sigma)
            grid[x,N] = fN(delta, np.array([grid[x+dd,N],grid[x-dd,N],\
				grid[x,N-dd]]),sigma)
            grid[0,x] = fN(delta, np.array([grid[0,x+dd],grid[0,x-dd],\
				grid[dd,x]]),sigma)
            grid[N,x] = fN(delta, np.array([grid[N,x+dd],grid[N,x-dd],\
				grid[N-dd,x]]),sigma)

            if wrap:
                grid[x,N] = grid[x,0]
                grid[N,x] = grid[0,x]

        x_vec = range(dd,N-dd+1,D)
        y_vec = range(D,N-dd+1,D)
	#interpolate offset interior grid points
        grid[dd:N-dd+1:D,D:N-dd+1:D] = [[fN(delta, np.array([grid[x,y+dd],\
			grid[x,y-dd],grid[x+dd,y],grid[x-dd,y]]), sigma)\
			for y in y_vec] for x in x_vec]

        x_vec = range(D,N-dd+1,D)
        y_vec = range(dd,N-dd+1,D)

        if x_vec != []:
            grid[D:N-dd+1:D,dd:N-dd+1:D] = [[fN(delta,\
					np.array([grid[x,y+dd], grid[x,y-dd],\
                    grid[x+dd,y], grid[x-dd,y]]), sigma)\
					for y in y_vec] for x in x_vec]

        if addition:          
            vec = range(0,N+1,D)
            grid[0:N+1:D,0:N+1:D] += [[delta*random.gauss(0.0,sigma)\
				for y in vec]for x in vec]

            vec = range(dd,N-dd+1,D)
            grid[dd:N-dd+1:D,dd:N-dd+1:D] += [[delta*random.gauss(0.0,sigma)\
				for y in vec]for x in vec]

        D=D/2
        dd=dd/2
    if(normalise):
        grid += np.amin(grid)*-1 + bounds[0]
        grid = (grid/np.amax(grid)) * bounds[1]
    return grid
