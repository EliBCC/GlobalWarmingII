# Linear version of Jupyter notebook for AI grader

# Import packages
import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import math
from matplotlib.animation import FuncAnimation
from IPython.display import HTML, display

# Grid and Variable Initialization
ncol = 3         # grid size (number of cells)
nrow = ncol

nSlices = 500         # maximum number of frames to show in the plot
ntAnim = 100          # number of time steps for each frame

horizontalWrap = True # determines whether the flow wraps around, connecting
                       # the left and right-hand sides of the grid, or whether
                       # there's a wall there. 
interpolateRotation = True
rotationScheme = "PlusMinus"   # "WithLatitude", "PlusMinus", "Uniform"
# rotationScheme = ""

# Note: the rotation rate gradient is more intense than the real world, so that
# the model can equilibrate quickly.

windScheme = ""  # "Curled", "Uniform"
initialPerturbation = "Tower"    # "Tower", "NSGradient", "EWGradient"
textOutput = False
plotOutput = True
arrowScale = 30

dT = 600 # seconds
# dT = 1000 # seconds
G = 9.8e-4 # m/s2, hacked (low-G) to make it run faster
HBackground = 4000 # meters

dX = 10.E3 # meters, small enough to respond quickly.  This is a very small ocean
# on a very small, low-G planet.  
dY = dX    # Grid is a square

dxDegrees = dX / 110.e3
flowConst = G  # 1/s2
dragConst = 1.E-6  # about 10 days decay time
meanLatitude = 30 # degrees

# Here's stuff you probably won't need to change

latitude = []
rotConst = []
windU = []
U = numpy.zeros((nrow, ncol+1))
V = numpy.zeros((nrow+1, ncol))
H = numpy.zeros((nrow, ncol+1))
dUdT = numpy.zeros((nrow, ncol))
dVdT = numpy.zeros((nrow, ncol))
dHdT = numpy.zeros((nrow, ncol))
dHdX = numpy.zeros((nrow, ncol+1))
dHdY = numpy.zeros((nrow, ncol))
dUdX = numpy.zeros((nrow, ncol))
dVdY = numpy.zeros((nrow, ncol))
rotV = numpy.zeros((nrow,ncol)) # interpolated to u locations
rotU = numpy.zeros((nrow,ncol)) #              to v

midCell = int(ncol/2)

# Set up reusable reset function to allow reruns of the script without restarting
def reset():
    latitude.clear()
    rotConst.clear()
    windU.clear()
    for irow in range(0,nrow):
        if rotationScheme == "WithLatitude":
            latitude.append( meanLatitude + (irow - nrow/2) * dxDegrees )
            rotConst.append( -7.e-5 * math.sin(math.radians(latitude[-1]))) # s-1
        elif rotationScheme == "PlusMinus":
            rotConst.append( -3.5e-5 * (1. - 0.8 * ( irow - (nrow-1)/2 ) / nrow )) # rot 50% +-
        elif rotationScheme == "Uniform":
            rotConst.append( -3.5e-5 ) 
        else:
            rotConst.append( 0 )

        if windScheme == "Curled":
            windU.append( 1e-8 * math.sin( (irow+0.5)/nrow * 2 * 3.14 ) ) 
        elif windScheme == "Uniform":
            windU.append( 1.e-8 )
        else:
            windU.append( 0 )
    global itGlobal
    itGlobal = 0

    U[:,:] = 0
    V[:,:] = 0
    H[:,:] = 0
    dUdT[:,:] = 0
    dVdT[:,:] = 0
    dHdT[:,:] = 0
    dHdX[:,:] = 0
    dHdY[:,:] = 0
    dUdX[:,:] = 0
    dVdY[:,:] = 0
    rotV[:,:] = 0
    rotU[:,:] = 0
        
    if initialPerturbation == "Tower":
        H[midCell,midCell] = 1
    elif initialPerturbation == "NSGradient":
        H[0:midCell,:] = 0.1
    elif initialPerturbation == "EWGradient":
        H[:,0:midCell] = 0.1

# Calculate dUdT based on formula provided
def longitudinalTimeDerivative(rotation, C_flow, dHdX_cell, C_drag, U_cell, C_wind):
    return rotation - C_flow * dHdX_cell - C_drag * U_cell + C_wind  

# Calculate dVdT based on formula provided
def latitudinalTimeDerivative(rotation, C_flow, dHdY_cell, C_drag, V_cell):
    return -rotation - C_flow * dHdY_cell - C_drag * V_cell

# Calculate dHdT based on formula provided
def verticalTimeDerivative(dUdX_cell, dVdY_cell, H_const, dX_const):
    return -(dUdX_cell + dVdY_cell) * H_const / dX_const

# Calculate interpolated rotation
def rotationInterpolation(U, V, rotU, rotV, rotConst):
    # Interpolate by averaging cells on either side of H
    U_interpolated = (U[:,1:] + U[:,:-1]) / 2
    V_interpolated = (V[1:,:] + V[:-1,:]) / 2

    # Multiply by rotation constant (varies only in the longitudinal direction)
    rotU_interpolated = numpy.array([rotConst]).T * U_interpolated
    rotV_interpolated = numpy.array([rotConst]).T * V_interpolated

    # Extend the interpolated rotations to allow for back-interpolation
    # Because we are using rotV to calculate U, and rotU to calculate V, 
    # we need to back-interpolate so that rotV is mapped on to the U positions, 
    # and vice versa
    if horizontalWrap:
        # If we have horizontal wrap, then the first column needs to be interpolated from the last column
        rotV_interpolated = numpy.hstack((rotV_interpolated[:, -1:], rotV_interpolated))
    else:
        # Otherwise, there is no preceding rotation, so use a value of 0
        rotV_interpolated = numpy.hstack((numpy.zeros((rotV_interpolated.shape[0], 1)), rotV_interpolated))

    # Longitudinally, there is never wrapping
    rotU_interpolated = numpy.vstack((numpy.zeros((1,rotU_interpolated.shape[1])), rotU_interpolated))

    # Back-interpolate
    rotU[:,:] = (rotU_interpolated[1:,:] + rotU_interpolated[:-1,:]) / 2
    rotV[:,:] = (rotV_interpolated[:,1:] + rotV_interpolated[:,:-1]) / 2


"""
This is the work-horse subroutine.  It steps forward in time, taking ntAnim steps of
duration dT.  
"""

def animStep():    

    global stepDump, itGlobal

    # Time Loop
    for it in range(0,ntAnim):
        # Here is where you need to build some code
    
        # Encode Longitudinal Derivatives Here
        # Calculate dHdX. We are mapping values from H to U, so we need to average the values of H on either side of U. In this case, those are H[i,j] and H[i,j-1]
        for i in range(len(dHdX)):
             for j in range(1, len(dHdX[0])):
                  dHdX[i, j] = (H[i, j] - H[i, j-1]) / dX
        if horizontalWrap:
             # We use a value of -2 here. -1 would be the ghost cells, which we don't want. We want the last column of real cells
             dHdX[:,0] = (H[:,0] - H[:,-2]) / dX
        else:
             dHdX[:,-1] = dHdX[:,0] = 0

        # Calculate dUdX. We are mapping U onto H, so we need to interpolate the values of U on either side of H. In this case, that's U[i,j+1] and U[i,j]
        for i in range(len(dUdX)):
             for j in range(len(dUdX[0])):
                  dUdX[i, j] = (U[i, j+1] - U[i, j]) / dX
        
        # Encode Latitudinal Derivatives Here
        # Calculate dHdY. We are mapping values from H to V, so we need to average the values of H on either side of V. In this case, those are H[i,j] and H[i-1,j]
        for i in range(1, len(dHdY)):
             for j in range(len(dHdY[0])):
                  dHdY[i, j] = (H[i, j] - H[i-1, j] ) / dY
     
        # Calculate dVdY. We are mapping V onto H, so we need to interpolate the values of V on either side of H. In this case, that's V[i+1,j] and V[i,j]
        for i in range(len(dVdY)):
              for j in range(len(dVdY[0])):
                   dVdY[i, j] = (V[i+1, j] - V[i, j]) / dY
        
        # Calculate the Rotational Terms Here
        
        if interpolateRotation:
             rotationInterpolation(U, V, rotU, rotV, rotConst)
        else:
          # Simple rotation calculation
          for i in range(len(rotU)):
               for j in range(len(rotU[0])):
                    rotU[i, j] = rotConst[i] * U[i, j]
                    rotV[i, j] = rotConst[i] * V[i, j]

        # Assemble the Time Derivatives Here
        for i in range(len(dUdT)):
               for j in range(len(dUdT[0])):
                    dUdT[i,j] = longitudinalTimeDerivative(rotV[i,j], flowConst, dHdX[i, j], dragConst, U[i, j], windU[i])
        
        for i in range(len(dVdT)):
            for j in range(len(dVdT[0])):
                    dVdT[i,j] = latitudinalTimeDerivative(rotU[i,j], flowConst, dHdY[i, j], dragConst, V[i, j])

        for i in range(len(dHdT)):
             for j in range(len(dHdT[0])):
                  dHdT[i,j] = verticalTimeDerivative(dUdX[i,j], dVdY[i,j], HBackground, dX)
        
        # Step Forward One Time Step
        # Apply time derivatives to value arrays
        for i in range(len(dUdT)):
             for j in range(len(dUdT[0])):
                  U[i, j] += dUdT[i, j] * dT 
        
        for i in range(len(dVdT)):
             for j in range(len(dVdT[0])):
                  V[i, j] += dVdT[i, j] * dT

        for i in range(len(dHdT)):
             for j in range(len(dHdT[0])):
                  H[i, j] += dHdT[i, j] * dT
    
        # Update the Boundary and Ghost Cells
        # Top and bottom flows are held at 0 as a boundary condition
        V[0,:] = V[-1,:] = 0

        if horizontalWrap:
             # Set the ghost cells to have the values from the left column, so that they are applied correctly in the next run
             U[:,-1] = U[:,0]
             H[:,-1] = H[:,0]
        else:
             # Hold the edges at 0 flow
             U[:,0] = U[:,-1] = 0

    #   Now you're done
    
    itGlobal = itGlobal + ntAnim

# Ensure all varaibles are fresh for this run
reset()
fig, ax = plt.subplots()
ax.set_title("H") 

hh = H[:, 0:ncol]
loc = tkr.IndexLocator(base=1, offset=1)
ax.xaxis.set_major_locator(loc)
ax.yaxis.set_major_locator(loc)
ax.grid(which='major', axis='both', linestyle='-')

hPlot = ax.imshow(hh, interpolation='nearest', clim=(-0.5, 0.5)) 

# Helper to compute arrow locations and vectors
def get_arrow_data():
    xx, yy, uu, vv = [], [], [], []
    # Horizontal arrows
    for irow in range(nrow):
        for icol in range(ncol):
            xx.append(icol - 0.5)
            yy.append(irow)
            uu.append(U[irow, icol] * arrowScale)
            vv.append(0)
    # Vertical arrows
    for irow in range(nrow):
        for icol in range(ncol):
            xx.append(icol)
            yy.append(irow - 0.5)
            uu.append(0)
            vv.append(-V[irow, icol] * arrowScale)
    return xx, yy, uu, vv

# Create initial quiver plot
xx, yy, uu, vv = get_arrow_data()
quiv = ax.quiver(xx, yy, uu, vv, color='white', scale=1)

def update(frame):
    # Step the simulation forward
    animStep()
    
    # Update height map (imshow)
    hh = H[:, 0:ncol]
    hPlot.set_array(hh)
    
    # Update vector fields (quiver)
    _, _, uu_new, vv_new = get_arrow_data()
    quiv.set_UVC(uu_new, vv_new)
    
    return hPlot, quiv

# 3. Build and display the animation
def init():
    return hPlot, quiv

# Using FuncAnimation to create a downloadable video, instead of a series of still frames
ani = FuncAnimation(
    fig,
    update,
    frames=nSlices,
    init_func=init,
    interval=100,
    blit=True
)

# Prevent duplicate static plot displaying in output
plt.close(fig)

# When this is moved to a Jupyter notebook, this will show an interactive video
if plotOutput is True:
    display(HTML(ani.to_jshtml()))

def textDump():
    print("time step ", itGlobal)    
    print("H", H)
    print("dHdX" )
    print( dHdX)
    print("dHdY" )
    print( dHdY)
    print("U" )
    print( U)
    print("dUdX" )
    print( dUdX)
    print("rotV" )
    print( rotV)
    print("V" )
    print( V)
    print("dVdY" )
    print( dVdY)
    print("rotU" )
    print( rotU)
    print("dHdT" )
    print( dHdT)
    print("dUdT" )
    print( dUdT)
    print("dVdT" )
    print( dVdT)

if textOutput is True:
    textDump()