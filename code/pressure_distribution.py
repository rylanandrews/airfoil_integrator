import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import random
import csv

import airfoil as af

class PressureDistribution:
    def __init__(self, alpha, airfoil, topPressures, bottomPressures, rhoWater, g=9.81):
        self.alpha = alpha
        self.airfoil: af.Airfoil = airfoil
        
        self.topPressures = topPressures.to_numpy()
        self.bottomPressures = bottomPressures.to_numpy()

        self.topAvg = np.full(9, np.nan)
        self.bottomAvg = np.full(9, np.nan)

        self.rhoWater = rhoWater
        self.g = g

    # Computes average pressure at a given station
    def getAvgPressure(self, stationNum, side):
        
        sum = 0
        count = 0

        if (side == 'top'):
            if (math.isnan(self.topAvg[stationNum])):
                for num in self.topPressures[:,stationNum]:
                    if (not math.isnan(num)):
                        sum += num
                        count += 1

                avg = sum / count
                # Convert to Pa
                avg = ((avg / 39.37) * self.g * self.rhoWater)
                self.topAvg[stationNum] = avg

            return self.topAvg[stationNum]
        
        if (side == 'bottom'):
            if (math.isnan(self.bottomAvg[stationNum])):
                for num in self.bottomPressures[:,stationNum]:
                    if (not math.isnan(num)):
                        sum += num
                        count += 1
                
                avg = sum / count
                # Convert to Pa
                avg = ((avg / 39.37) * self.g * self.rhoWater)
                self.bottomAvg[stationNum] = avg

            return self.bottomAvg[stationNum]

        return -1
            
    
    # Returns standard deviation in inH2O
    def getStdDevPressure(self, stationNum, side):
        # Unconvert from Pa
        avg = self.getAvgPressure(stationNum, side) * (39.37 / (self.g * self.rhoWater))

        sum = 0
        count = 0

        if (side == 'top'):
            for num in self.topPressures[:,stationNum]:
                if (not math.isnan(num)):
                    sum += (num - avg)**2
                    count += 1
        
        if (side == 'bottom'):
            for num in self.bottomPressures[:,stationNum]:
                if (not math.isnan(num)):
                    sum += (num - avg)**2
                    count += 1

        return math.sqrt( (1/(count-1)) * sum )
    
    def getPressureUncertainty(self, stationNum, side, verbose = False):
        U_rhoWater = 11
        U_g = 0.01

        U_sys = 0.04

        avg = self.getAvgPressure(stationNum, side) * (39.37 / (self.g * self.rhoWater))
        stdDev = self.getStdDevPressure(stationNum, side)

        U_rand = 1.96 * stdDev
        U_h2o = math.sqrt(U_sys**2 + U_rand**2)

        term1 = U_rhoWater * self.g * (avg / 39.37)
        term2 = U_g * self.rhoWater * (avg / 39.37)
        term3 = (U_h2o / 39.37) * self.rhoWater * self.g

        if verbose:
            print("Random: " + str(U_rand))
            print("Systematic: " + str(U_sys))
            print("Total: " + str(U_h2o))
            print("Density Term: " + str(term1))
            print("Gravity Term 2: " + str(term2))
            print("Pressure Term 3: " + str(term3))


        return math.sqrt(term1**2 + term2**2 + term3**2)

    
    # Plots the pressure distribution using the normalal vectors to the airfoil
    def plotPressureDist(self):

        topnormals = self.airfoil.getPortNorms(side='top')
        for i in range(0,len(topnormals)):
            avgPressure = -self.getAvgPressure(stationNum=i, side='top')
            topnormals[i,:] = avgPressure * topnormals[i,:]

        bottomnormals = self.airfoil.getPortNorms(side='bottom')
        for i in range(0,len(bottomnormals)):
            avgPressure = -self.getAvgPressure(stationNum=i, side='bottom')
            bottomnormals[i,:] = avgPressure * bottomnormals[i,:]


        plt.quiver(self.airfoil.xPorts, self.airfoil.yPorts, topnormals[:,0], topnormals[:,1])
        plt.quiver(self.airfoil.xPorts, -self.airfoil.yPorts, bottomnormals[:,0], bottomnormals[:,1])
        plt.plot(self.airfoil.xCoords, self.airfoil.yCoords)
        plt.plot(self.airfoil.xPorts, self.airfoil.yPorts, '*')
        plt.plot(self.airfoil.xPorts, -self.airfoil.yPorts, '*')

        plt.title("Angle of Attack (deg): " + str(self.alpha))
        plt.axis('equal')
        plt.show()

    def sumForcesAirfoil(self, N=80):
        pTop = []
        pBottom = []

        for i in range(9):
            pTop.append(self.getAvgPressure(i, 'top'))
            pBottom.append(self.getAvgPressure(i, 'bottom')) 

        return self.sumForcesAirfoilHelper(pTop, pBottom, N)

    # Sum x and y components of force
    def sumForcesAirfoilHelper(self, pTop, pBottom, N=80):
        sumX = 0
        sumY = 0

        x = np.linspace(0, 1, num=N)

        currentPort = 0

        for i in range(len(x)-1):
            x1 = x[i]
            x3 = x[i+1]
            x2 = (x3 + x1) / 2
            
            if (currentPort < 8):
                if ( abs(self.airfoil.xPorts[currentPort] - x2) > abs(self.airfoil.xPorts[currentPort+1] - x2) ):
                    currentPort += 1

            forceTop = self.computeForce(x1, x2, x3, 'top', pTop[currentPort])
            forceBottom = self.computeForce(x1, x2, x3, 'bottom', pBottom[currentPort])

            sumX = sumX + forceTop[0] + forceBottom[0]
            sumY = sumY + forceTop[1] + forceBottom[1]

        return [-sumX, -sumY]
    
    def computeForce(self, x1, x2, x3, side, pressure):
        y1 = self.airfoil.Y(x=x1)
        y3 = self.airfoil.Y(x=x3)

        area = math.dist([x1 * self.airfoil.chordLength, y1 * self.airfoil.chordLength], [x3 * self.airfoil.chordLength, y3 * self.airfoil.chordLength])

        normal = self.airfoil.getNorm(x=x2, side=side)

        return normal * pressure * area
    
    def sumForcesLD(self, N=80):
        return AC2VC(self.sumForcesAirfoil(N), self.alpha)
    
    def forceUncertainty(self, numIter=4000):
        U_alpha = 0.5

        # Find averages and uncertainties
        pTop = []
        pBottom = []
        U_pTop = []
        U_pBottom = []

        for i in range(9):
            pTop.append(self.getAvgPressure(i, 'top'))
            pBottom.append(self.getAvgPressure(i, 'bottom'))
            U_pTop.append(self.getPressureUncertainty(i, 'top'))
            U_pBottom.append(self.getPressureUncertainty(i, 'bottom'))
         
        
        # Repeat for simulation num:
        #   Generate new pressures based on uncertainty
        #   Sum forces using new pressures
        #   Generate new angle of attack based on uncertainty
        #   Compute lift and drag from rotation
        #   Save lift, drag to a list

        liftSims = []
        dragSims = []

        for i in range(numIter):
            pTopTemp = []
            pBottomTemp = []

            for i in range(9):
                pTopTemp.append(norm.ppf(random.random(), loc=pTop[i], scale=U_pTop[i]))
                pBottomTemp.append(norm.ppf(random.random(), loc=pBottom[i], scale=U_pBottom[i]))

            vec = self.sumForcesAirfoilHelper(pTopTemp, pBottomTemp)

            alphaTemp = norm.ppf(random.random(), loc=self.alpha, scale=U_alpha)

            [drag, lift] = AC2VC(vec, alphaTemp)

            liftSims.append(lift)
            dragSims.append(drag)

        # Compute std. deviation in lift and drag
        sumLift = 0
        sumDrag = 0
        for i in range(numIter):
            sumLift += liftSims[i]
            sumDrag += dragSims[i]

        avgLift = sumLift / numIter
        avgDrag = sumDrag / numIter

        sdSumLift = 0
        sdSumDrag = 0
        for i in range(numIter):
            sdSumLift += (liftSims[i] - avgLift)**2
            sdSumDrag += (dragSims[i] - avgDrag)**2

        stdDevLift = math.sqrt( (1/(numIter-1)) * sdSumLift )
        stdDevDrag = math.sqrt( (1/(numIter-1)) * sdSumDrag )

        with open('lift_sims.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(liftSims)

        with open('drag_sims.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(dragSims)
        
        return [1.96 * stdDevLift, 1.96 * stdDevDrag]
    
    def computeMoment(self, xm, N=80):
        pTop = []
        pBottom = []

        for i in range(9):
            pTop.append(self.getAvgPressure(i, 'top'))
            pBottom.append(self.getAvgPressure(i, 'bottom'))

        return self.computeMomentHelper(xm, pTop, pBottom, self.alpha, N)
    
    def computeMomentHelper(self, xm, pTop, pBottom, alpha, N=80):
        totMoment = 0

        x = np.linspace(0, 1, num=N)

        currentPort = 0

        for i in range(len(x)-1):
            x1 = x[i]
            x3 = x[i+1]
            x2 = (x3 + x1) / 2

            y2 = self.airfoil.Y(x2)
            
            if (currentPort < 8):
                if ( abs(self.airfoil.xPorts[currentPort] - x2) > abs(self.airfoil.xPorts[currentPort+1] - x2) ):
                    currentPort += 1

            forceTop = self.computeForce(x1, x2, x3, 'top', pTop[currentPort])
            forceBottom = self.computeForce(x1, x2, x3, 'bottom', pBottom[currentPort])

            forceTopRot = AC2VC(forceTop, alpha)
            forceBottomRot = AC2VC(forceBottom, alpha)

            totMoment = totMoment + (forceTopRot[1] * self.airfoil.chordLength * (x2-xm)) + (forceBottomRot[1] * self.airfoil.chordLength * (x2-xm))

        return totMoment
    
    def momentUncertainty(self, xm, numIter=4000):
        U_alpha = 0.5

        # Find averages and uncertainties
        pTop = []
        pBottom = []
        U_pTop = []
        U_pBottom = []

        for i in range(9):
            pTop.append(self.getAvgPressure(i, 'top'))
            pBottom.append(self.getAvgPressure(i, 'bottom'))
            U_pTop.append(self.getPressureUncertainty(i, 'top'))
            U_pBottom.append(self.getPressureUncertainty(i, 'bottom'))
         
        # Repeat for simulation num:
        #   Generate new pressures based on uncertainty
        #   Sum forces using new pressures
        #   Generate new angle of attack based on uncertainty
        #   Compute lift and drag from rotation
        #   Save lift, drag to a list

        momentSims = []

        for i in range(numIter):
            pTopTemp = []
            pBottomTemp = []

            for i in range(9):
                pTopTemp.append(norm.ppf(random.random(), loc=pTop[i], scale=U_pTop[i]))
                pBottomTemp.append(norm.ppf(random.random(), loc=pBottom[i], scale=U_pBottom[i]))

            alphaTemp = norm.ppf(random.random(), loc=self.alpha, scale=U_alpha)

            momentTemp = self.computeMomentHelper(xm, pTopTemp, pBottomTemp, alphaTemp)

            momentSims.append(momentTemp)

        # Compute std. deviation in lift and drag
        sum = 0
        for i in range(numIter):
            sum += momentSims[i]

        avg = sum / numIter

        sdSum = 0
        for i in range(numIter):
            sdSum += (momentSims[i] - avg)**2

        stdDevMom = math.sqrt( (1/(numIter-1)) * sdSum )
        
        return 1.96 * stdDevMom

# Rotates a vector from one basis to another by angle alpha
def AC2VC(vector, alpha):
        x = vector[0]
        y = vector[1]
        angle = math.radians(alpha)

        print("Before rotation: [" + str(x) + ", " + str(y) + "]")

        rot = [x*math.cos(angle) + y*math.sin(angle), -x*math.sin(angle) + y*math.cos(angle)]

        print("After rotation: [" + str(rot[0]) + ", " + str(rot[1]) + "]")

        
        return rot








        