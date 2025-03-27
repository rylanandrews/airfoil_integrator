# Code for analysis of airfoils in Lab 5
# Author: Rylan Andrews
# Date: 2025-03-04

import pandas as pd
import csv
import time

import airfoil as af
import pressure_distribution as pdist


# Constants for convenience
N = 80
xm = 0.305
numIter = 4000

def main():
    # Read in airfoil data
    airfoilCoords = pd.read_csv("code\\naca0012.csv", header=None)
    portCoords = pd.read_csv("code\\port_chord_pos.csv", header=None)

    # Generate the airfoil
    airfoil = af.Airfoil(chordLength=(4/39.37), airfoilCoords=airfoilCoords, portCoords=portCoords)

    liftDragVectors = []
    ld_uncertainties = []
    moments = []
    mom_uncertainties = []

    # Load data
    pressureData = pd.read_csv("combined_data\\merged_data_a0.csv", header=None)
    pressureDistributions = {
        '0': pdist.PressureDistribution(alpha=0, airfoil=airfoil, topPressures=pressureData.iloc[:,0:9], bottomPressures=pressureData.iloc[:,0:9], rhoWater=1000)
    }
    liftDragVectors.append(pressureDistributions['0'].sumForcesLD(N = N))
    #ld_uncertainties.append(pressureDistributions['0'].forceUncertainty(numIter))
    moments.append(pressureDistributions['0'].computeMoment(xm=xm))
    #mom_uncertainties.append(pressureDistributions['0'].momentUncertainty(xm, numIter))
    #pressureDistributions['0'].plotPressureDist()

    

    for i in [2, 4, 6, 8, 9, 10, 11, 12, 13, 14]:
        pressureData = pd.read_csv("combined_data\\merged_data_a" + str(i) + ".csv")
        pressureDistributions[str(i)] = pdist.PressureDistribution(alpha=i, airfoil=airfoil, topPressures=pressureData.iloc[:,0:9], bottomPressures=pressureData.iloc[:,9:18], rhoWater=1000)
        liftDragVectors.append(pressureDistributions[str(i)].sumForcesLD(N))
        #ld_uncertainties.append(pressureDistributions[str(i)].forceUncertainty(numIter))
        moments.append(pressureDistributions[str(i)].computeMoment(xm=xm))
        #mom_uncertainties.append(pressureDistributions[str(i)].momentUncertainty(xm, numIter))
        #pressureDistributions[str(i)].plotPressureDist()

    
    with open('lift_drag.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(liftDragVectors)

    with open('moments.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(moments)

    with open('ld_uncertainty.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(ld_uncertainties)

    with open('mom_uncertainties.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(mom_uncertainties)

    print("Success!")

startTime = time.time()

main()

endTime = time.time()

print("Elapsed time (min): " + str((endTime - startTime)/60))
