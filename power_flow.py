#busdata must be in order in excel file, buses must start at 1 and go straight up to max


import numpy as np
import pandas as pd

class PowerFlow():

    def __init__(self, maxIterations, tolerance):

        self.maxIterations = maxIterations
        self.currIterations = 0
        
        self.tolerance = tolerance

        self.BusData = 0
        self.LineData = 0
        self.numBusses = 0
        self.numConnections = 0
        
        self.admittanceReal = 0
        self.admittanceImag = 0
        
        #self.equationMatrix = 0
        self.PEquationMatrix = 0
        self.QEquationMatrix = 0
        self.implicitEquationMatrix = 0
        
        self.deltaMatrix = 0
        self.inverseJacobian = 0

        self.voltages = 0
        self.angles = 0
        self.busType = 0

        self.convergenceHistory = 0

    #pulls data from excel file into class
    def readFromFile(self, fileName):
        self.BusData = pd.read_excel(fileName, sheet_name = 0)
        self.LineData = pd.read_excel(fileName, sheet_name = 1)

        self.numBusses = len(self.BusData) # counts the number of buses
        self.numConnections = len(self.LineData) # counts the number of connections

        #build matrices that will not change over NR iterations
        self.buildAdmittanceMatrix()
        self.getVoltages()
        self.getBusType()

        #update equations matrix with processed values to prepare for 1st NR iteration
        self.updateEquationMatrix()

    def buildAdmittanceMatrix(self):
        admittanceComplex = np.zeros((self.numBusses, self.numBusses), dtype=np.complex_)
        self.admittanceReal = np.zeros((self.numBusses, self.numBusses))
        self.admittanceImag = np.zeros((self.numBusses, self.numBusses))

        #iterate through connection, and for each add to every cell of admittance matrix it affects
        for i in range(self.numConnections):
            node1 = self.LineData['From'][i]
            node2 = self.LineData['To'][i]

            #add the admittance between nodes to corresponding off diagonal rows
            #then add the admittance between nodes and 1/2 the susceptance to the diagonal rows

            resistance1to2 = self.LineData['Rtotal, p.u.'][i]
            reactance1to2 = self.LineData['Xtotal, p.u.'][i]
            reactanceShunt = self.LineData['Btotal, p.u.'][i]

            admittance1to2 = 1/(resistance1to2 + 1j * reactance1to2)
            susceptance1and2 = (1j * reactanceShunt) / 2

            admittanceComplex[node1 - 1][node1 - 1] += ( admittance1to2 + susceptance1and2 )
            admittanceComplex[node2 - 1][node2 - 1] += ( admittance1to2 + susceptance1and2 )

            admittanceComplex[node1 - 1][node2 - 1] -= ( admittance1to2 )
            admittanceComplex[node2 - 1][node1 - 1] -= ( admittance1to2 )

        self.admittanceReal = np.real(admittanceComplex)
        self.admittanceImag = np.imag(admittanceComplex)
        
    def getVoltages(self): #voltages AND angles
        self.voltages = self.BusData['V Set']
        self.angles = np.zeros_like(self.voltages)

    def getBusType(self):
        self.busType = self.BusData['Type']

    def newtonRaphsonIteration(self):
        #update all these values
        self.currIterations += 1
        
        self.updateinverseJacobian()
        self.updateDeltaMatrix()
        self.updateVoltages()
        self.updateEquationMatrix()

        return self.deltaMatrix, self.equationMatrix

    def updateInverseJacobian(self):
        numPQ = self.busType('D')
        numPV = self.busType('G')
        jacobian = np.zeros(self.numBusses + numPQ, self.numBusses + numPQ)
        jacobianSum = 0

        # H matrix
        for i in range(1, self.numBusses):
            for k in range(1, self.numBusses):
                if i != k:
                    jacobian[k - 1][i - 1] = self.voltage[k] * self.voltage[i] * (self.admittanceReal[k][i] * np.sin(
                                             self.angle[k] - self.angle[i]) - self.admittanceImag[k][i] * np.cos(
                                             self.angle[k] - self.angle[i]))
                
                    jacobian[i - 1][i - 1] -= jacobian[k][i]

        # M matrix
        for i in range(numPV + 1, self.numBusses):
            for k in range(1, self.numBusses):
                if i != k:
                    jacobian[k - 1][i + numPQ - 1] = self.voltage[k] * (self.admittanceReal[k] * np.cos(
                                                      self.angle[k] - self.angle[i]) + self.admittanceImag * np.sin(
                                                      self.angle[k] - self.angle[i]))
                
                    jacobian[i + numPQ - 1][i + numPQ - 1] += jacobian[k - 1][i + numPQ] - 1 +  2 * self.admittanceReal[k][k] * self.voltages[k]

        # N matrix
        for i in range(1, self.numBusses):
            for k in range(numPV + 1, self.numBusses):
                if i != k:
                    jacobian[k + numPQ - 1][i - 1] = self.voltage[k] * self.voltage[i] * (np.negative(self.admittanceReal[k][i]) * np.cos(
                                                 self.angle[k] - self.angle[i]) - self.admittanceImag[k][i] * np.sin(
                                                 self.angle[k] - self.angle[i]))
                
                    jacobian[i - 1][i - 1] -= jacobian[k + numPQ - 1][i - 1]

        # L matrix
        for i in range(numPV + 1, self.numBusses):
            for k in range(numPV + 1, self.numBusses):
                if i != k:
                    jacobian[k + numPQ - 1][i + numPQ - 1] = self.voltage[k] * (self.admittanceReal[k][i] * np.sin(
                                                     self.angle[k] - self.angle[i]) - self.admittanceImag[k][i] * np.cos(
                                                     self.angle[k] - self.angle[i]))

                    jacobianSum += jacobian[k + numPQ - 1][i + numPQ - 1] 
                    
            jacobian[i + numPQ - 1][i + numPQ - 1] = np.negative(2) * self.admittanceImag[k][k] * self.voltages[k] + jacobianSum

        
        self.inverseJacobian = np.linalg.inv(jacobian)

    def updateDeltaMatrix(self):
        return
        #self.deltaMatrix = np.dot( (-1 * self.inverseJacobian), self.equationMatrix)

    def updateVoltages(self):
        return
        #these won't be the same shape, I don't think. Will not be this simple
        #self.voltages += self.deltaMatrix

    def updateEquationMatrix(self):
        self.PEquationList = np.zeros(self.numBusses)
        self.QEquationList = np.zeros(self.numBusses)
        #self.implicitEquationMatrix = 

        for i in range(self.numBusses):
            for j in range(self.numBusses):
                self.PEquationList[i][j] += ( self.voltages[i] * self.voltages[j] * ( self.admittanceReal[i][j]*np.cos(self.angles[i] - self.angles[j]) +self.admittanceImag[i][j]*np.sin(self.angles[i] - self.angles[j] ) ) )

                self.QEquationList[i][j] += ( self.voltages[i] * self.voltages[j] * ( self.admittanceReal[i][j]*np.sin(self.angles[i] - self.angles[j]) -self.admittanceImag[i][j]*np.cos(self.angles[i] - self.angles[j] ) ) )

        print(self.EquationMatrix)
        
        #for loop through matrix and update each equation according to power flow slides
        #if the equation has a known P or Q, then also add it to the implicit equations list
        #will have to look at the type of bus to do this, I think? PQ buses produce implicit equations?
        #also seems like we will end up with too many implicit equations, how do we pick which ones to use?
