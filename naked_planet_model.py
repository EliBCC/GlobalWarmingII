"""
The goal is to numerically simulate how the planetary temperature of a naked planet would change through 
time as it approaches equilibrium (the state at which it stops changing, which we calculated before).  
The planet starts with some initial temperature.  The “heat capacity” (units of Joules / m2 K) of the planet 
is set by a layer of water which absorbs heat and changes its temperature.  If the layer is very thick, it 
takes a lot more heat (Joules) to change the temperature. The differential equation you are going to solve is 

dHeatContent/dt = L*(1-alpha)/4 - epsilon * sigma * T^4

where the heat content is related to the temperature by the heat capacity   

T[K] = HeatContent [J/m2] / HeatCapacity [J/m2 K]

The numerical method is to take time steps, extrapolating the heat content from one step to the next using the 
incoming and outgoing heat fluxes, same as you would balance a bank account by adding all the income and 
subtracting all the expenditures over some time interval like a month.  The heat content "jumps" from the 
value at the beginning of the time step, to the value at the end, by following the equation

HeatContent(t+1) = HeatContent(t) + dHeatContent/dT * TimeStep

This scheme only works if the time step is short enough that nothing too huge happens through the course of the 
time step.  

Set the model up with the following constants:

timeStep = 100           # years
waterDepth = 4000        # meters
L = 1350                 # Watts/m2
albedo = 0.3
epsilon = 1
sigma = 5.67E-8          # W/m2 K4

"""

import numpy as np
import matplotlib.pyplot as plt

steps = 500
timeStep = 10           # years
waterDepth = 4000        # meters
L = 1350                 # Watts/m2
albedo = 0.3
epsilon = 1
sigma = 5.67E-8          # W/m2 K4
initialTemp = 0          # K

heatCapacityCubicMWater = 4184000   # J/m^3
heatCapacityEarth = waterDepth * heatCapacityCubicMWater  # J/m^2

secondsPerYear = 60 * 60 * 24 * 365 # s/year

def temperature(heat_content):
    return heat_content / heatCapacityEarth

def incomingSolar():
    return L * (1 - albedo) / 4

def outgoingIR(temp):
    return epsilon * sigma * pow(temp, 4)

def totalHeatFluxWatts(temp):
    return incomingSolar() - outgoingIR(temp)

def totalHeatFluxJoules(temp):
    return totalHeatFluxWatts(temp) * secondsPerYear * timeStep

def main():
    # Initialize arrays
    times = np.arange(steps) * timeStep
    temps = np.empty(steps)
    heat_content = 0


    # Set first row (calculate initial heat content from initial temp)
    temps[0] = initialTemp
    heat_content = temps[0] * heatCapacityEarth

    # Simulate
    for s in range(1, steps):
        heat_content += totalHeatFluxJoules(temps[s - 1])
        temps[s] = temperature(heat_content)

    # Plot table
    plt.figure(figsize=(8,5))
    plt.plot(times, temps, marker='o', linestyle='-', color='b', label='Temperature vs Time')

    plt.title('Naked Planet Model')
    plt.xlabel('Time (years)')
    plt.ylabel('Temperature (K)')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()

    plt.savefig('naked_planet_model.png', dpi=300)
    print("Plot successfully saved!")

if __name__ == "__main__":
    steps = int(input(""))
    main()