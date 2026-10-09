# Simpified script for autograder

import numpy
import math

# Grid and Variable Initialization
ncol = 3         # grid size (number of cells)
nrow = ncol
ncol = 3
nSlices, iRowOut, iColOut = input("").split()
nSlices, iRowOut, iColOut = [ int(nSlices), int(iRowOut), int(iColOut) ]         
ntAnim = 1          
horizontalWrap = True 
interpolateRotation = True  # or True, either works
textOutput = False
plotOutput = False

dT = 600    # seconds
G = 9.8e-4  # artificially low to allow a long time step
rotationScheme = "PlusMinus"

# Note: the rotation rate gradient is more intense than the real world, so that
# the model can equilibrate quickly.

windScheme = ""  # "Curled", "Uniform"
initialPerturbation = "Tower"    # "Tower", "NSGradient", "EWGradient"
textOutput = False
plotOutput = True
arrowScale = 30

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


def longitudinalTimeDerivative(rotation, C_flow, dHdX_cell, C_drag, U_cell, C_wind):
    return rotation - C_flow * dHdX_cell - C_drag * U_cell + C_wind  

def latitudinalTimeDerivative(rotation, C_flow, dHdY_cell, C_drag, V_cell):
    return -rotation - C_flow * dHdY_cell - C_drag * V_cell

def verticalTimeDerivative(dUdX_cell, dVdY_cell, H_const, dX_const):
    return -(dUdX_cell + dVdY_cell) * H_const / dX_const

def rotationInterpolation(U, V, rotU, rotV, rotConst):
    U_interpolated = (U[:,1:] + U[:,:-1]) / 2
    V_interpolated = (V[1:,:] + V[:-1,:]) / 2

    rotU_interpolated = numpy.array([rotConst]).T * U_interpolated
    rotV_interpolated = numpy.array([rotConst]).T * V_interpolated

    if horizontalWrap:
        rotV_interpolated = numpy.hstack((rotV_interpolated[:, -1:], rotV_interpolated))
    else:
        rotV_interpolated = numpy.hstack((numpy.zeros((rotV_interpolated.shape[0], 1)), rotV_interpolated))

    rotU_interpolated = numpy.vstack((numpy.zeros((1,rotU_interpolated.shape[1])), rotU_interpolated))

    rotU[:,:] = (rotU_interpolated[1:,:] + rotU_interpolated[:-1,:]) / 2
    rotV[:,:] = (rotV_interpolated[:,1:] + rotV_interpolated[:,:-1]) / 2

"""
This is the work-horse subroutine.  It steps forward in time, taking ntAnim steps of
duration dT.  
"""

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
        for i in range(len(dHdX)):
             for j in range(1, len(dHdX[0])):
                  dHdX[i, j] = (H[i, j] - H[i, j-1]) / dX
        if horizontalWrap:
             dHdX[:,0] = (H[:,0] - H[:,-2]) / dX
        else:
             dHdX[:,-1] = dHdX[:,0] = 0

        for i in range(len(dUdX)):
             for j in range(len(dUdX[0])):
                  dUdX[i, j] = (U[i, j+1] - U[i, j]) / dX
        
        # Encode Latitudinal Derivatives Here
        for i in range(1, len(dHdY)):
             for j in range(len(dHdY[0])):
                  dHdY[i, j] = (H[i, j] - H[i-1, j] ) / dY
     

        for i in range(len(dVdY)):
              for j in range(len(dVdY[0])):
                   dVdY[i, j] = (V[i+1, j] - V[i, j]) / dY
        
        # Calculate the Rotational Terms Here
        
        if interpolateRotation:
             rotationInterpolation(U, V, rotU, rotV, rotConst)
        else:
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
        V[0,:] = V[-1,:] = 0
        if horizontalWrap:
             U[:,-1] = U[:,0]
             H[:,-1] = H[:,0]
        else:
             U[:,0] = U[:,-1] = 0

    #   Now you're done
    
    itGlobal = itGlobal + ntAnim



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

for _ in range(nSlices):
    animStep()
     
print(H[iRowOut,iColOut],dHdT[iRowOut,iColOut],U[iRowOut,iColOut],V[iRowOut,iColOut],rotU[iRowOut,iColOut])
