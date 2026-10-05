# Simplified Python script version of the assignment
# Required for autograder
import numpy as np

nX = 10                # number of grid points
domainWidth = 1e6      # meters
timeStep = 100         # years
flowParam = 1e4        # m horizontal / yr
snowFall = 0.5         # m / y
dX = domainWidth / nX  # m / grid point

def calcFlows(flows):
    for i in range (nX + 2 - 1):
        flows[i] = ( elevations[i] - elevations[i+1] ) / dX * flowParam  * \
            ( elevations[i]+elevations[i+1] ) / 2 / dX

def calcElevations(flows, elevations):
    for i in range (1, nX + 2 - 1):
        elevations[i] += ( snowFall + flows[i-1] - flows[i] ) * timeStep

if __name__ == "__main__":
    elevations = np.zeros(nX+2)
    flows = np.zeros(nX+1)
    nYears = float( input('') )
    for _ in np.arange(0, nYears, timeStep):
        calcFlows(flows)
        calcElevations(flows, elevations)
    print(elevations[int(len(elevations) / 2)])