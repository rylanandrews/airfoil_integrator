import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class Airfoil:
    def __init__(self, chordLength, airfoilCoords, portCoords):
        # Airfoil chord length (in meters)
        self.chordLength = chordLength

        # Given coordinates of airfoil surfaces
        self.xCoords = airfoilCoords.iloc[:,0].to_numpy()
        self.yCoords = airfoilCoords.iloc[:,1].to_numpy()

        # Given coordinates of pressure ports
        self.xPorts = portCoords.iloc[:,0].to_numpy()
        self.yPorts = self.getYPorts()

    # Returns a y position for a given x, even if it's not one of the nominal coordinates
    def Y(self, x, negative = False):
        x1 = 0
        x2 = 1
        y1 = 0
        y2 = 0

        i1 = 0
        i2 = 0

        for i in range(len(self.xCoords)):
            if (self.xCoords[i] < x and self.xCoords[i] > x1):
                x1 = self.xCoords[i]
                y1 = self.yCoords[i]
                i1 = i
            
            if (self.xCoords[i] > x and self.xCoords[i] < x2):
                x2 = self.xCoords[i]
                y2 = self.yCoords[i]
                i2 = 0
        
        if (x1 == x2):
            x2 = self.xCoords[i2 + 1]
            y2 = self.yCoords[i2 + 1]
            
        y = y1 + ((y2 - y1) / (x2 - x1)) * (x - x1)

        if negative:
            return -y
        
        return y

    # Plots the nominal coordinates
    def plot(self):
        plt.plot(self.xCoords, self.yCoords, 'o')
        plt.axis('equal')
        plt.show()

    # Plots N interpolated points
    def plotN(self, N):
        xN = np.linspace(0, 1, num=N)
        yN = np.empty(N)

        for i in range(N):
            yN[i] = self.Y(xN[i])

        plt.plot(xN, yN, 'o')
        plt.plot(self.xCoords, self.yCoords, '*')
        plt.axis('equal')
        plt.show()

    # Computes the Y positions of the pressure ports
    def getYPorts(self):
        yPorts = np.empty((len(self.xPorts), 1))

        for i in range(len(self.xPorts)):
            yPorts[i] = self.Y(self.xPorts[i])

        return yPorts

    # Plots the locations of the pressure ports
    def plotPorts(self):
        yPorts = self.getYPorts()

        plt.plot(self.xCoords, self.yCoords)
        plt.plot(self.xPorts, yPorts, '*')
        plt.axis('equal')
        plt.show()

    # Saves the x and y locations of the pressure ports
    def savePorts(self, fileName):
        yPorts = self.getYPorts()
        portXY = np.concatenate((self.xPorts, yPorts))
        print(portXY)
        df = pd.DataFrame(portXY)
        df.to_csv(fileName)

    # Returns [x, y] normal vector to the surface
    def getNorm(self, x, side):
        # Find slope at x, y
        x1 = x - 0.001
        x2 = x + 0.001
        if x1 < 0:
            x1 = 0
        if x2 > 1:
            x2 = 1

        y1 = self.Y(x1)
        y2 = self.Y(x2)

        m = (y2 - y1) / (x2 - x1)

        # Compute vector
        v = np.array([1, -1/m])

        # Make sure vector points out from surface (positive y)
        if (v[1] < 0):
            v = -v

        if (side == 'bottom'):
            v[1] = -v[1]

        # Normalize
        v_mag = np.linalg.norm(v)
        if v_mag == 0:
            return v
        return v / v_mag
    
    def getPortNorms(self, side):
        
        norms = np.empty((0,2))

        for i in range(len(self.xPorts)):
            newNorm = self.getNorm(self.xPorts[i], side)
            norms = np.concatenate((norms, newNorm.reshape(1,-1)), axis=0)

        return norms

    # Plots the normal vectors at each port
    def plotPortNorms(self):
        norms = self.getPortNorms(side='top')

        plt.quiver(self.xPorts, self.yPorts, norms[:,0], norms[:,1])
        plt.plot(self.xCoords, self.yCoords)
        plt.plot(self.xPorts, self.yPorts, '*')
        plt.axis('equal')
        plt.show()


