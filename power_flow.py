#busdata must be in order in excel file, buses must start at 1 and go straight up to max

#equations list is correct for notes example
#jacobian is not correct for notest example

#i don't think anything is correct for basecase

#need to add in power from generator for those nodes
#check over how power is gathered... not doing it correctly in the equations


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
        
        self.PEquationList = 0
        self.QEquationList = 0
        self.implicitEquationList = 0
        
        self.deltaMatrix = 0
        self.inverseJacobian = 0

        self.voltages = 0
        self.angles = 0
        self.busType = 0

        self.busMap = 0

        self.convergenceHistory = 0

    #pulls data from excel file into class
    def readFromFile(self, fileName):
        self.BusData = pd.read_excel(fileName, sheet_name = 0)
        self.BusData = self.BusData.sort_values( by='Type', ascending = False) #sort data by bus type

        self.LineData = pd.read_excel(fileName, sheet_name = 1) #need to sort lineData?

        self.numBusses = len(self.BusData) # counts the number of buses
        self.numConnections = len(self.LineData) # counts the number of connections

        #make dictionary mapping new, sorted values to old ones
        #need to make a new bus number system with buses in order,
            #but have to remember the old system and bring that back in when needed
            #also need to translate from line data, that data will refer to old line numbers. Reorder line data?

        self.busMap = np.array(self.BusData.index)

        #build the matrices that will not change over NR iterations
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

            #these lines pull the from and to nodes listen in the file
            node1 = self.LineData['From'][i]
            node2 = self.LineData['To'][i]

            #these indexes store where the admittance for the specific node should go in the matrix
            node1Index = np.where(self.busMap == (node1 - 1))[0][0] - 1
            node2Index = np.where(self.busMap == (node2 - 1))[0][0] - 1

            #add the admittance between nodes to corresponding off diagonal rows
            #then add the admittance between nodes and 1/2 the susceptance to the diagonal rows

            resistance1to2 = self.LineData['Rtotal, p.u.'][i]
            reactance1to2 = self.LineData['Xtotal, p.u.'][i]
            reactanceShunt = self.LineData['Btotal, p.u.'][i]

            admittance1to2 = 1/(resistance1to2 + 1j * reactance1to2)
            susceptance1and2 = (1j * reactanceShunt) / 2

            admittanceComplex[node1Index][node1Index] += ( admittance1to2 + susceptance1and2 )
            admittanceComplex[node2Index][node2Index] += ( admittance1to2 + susceptance1and2 )

            admittanceComplex[node1Index][node2Index] -= ( admittance1to2 )
            admittanceComplex[node2Index][node1Index] -= ( admittance1to2 )

        self.admittanceReal = np.real(admittanceComplex)
        self.admittanceImag = np.imag(admittanceComplex)
        
    def getVoltages(self): #voltages AND angles
        self.voltages = np.array(self.BusData['V Set'])
        self.angles = np.zeros_like(self.voltages)

    def getBusType(self):
        self.busType = np.array(self.BusData['Type'])

    def newtonRaphsonIteration(self):
        #update all these values
        self.currIterations += 1
        
        self.updateInverseJacobian()
        self.updateDeltaMatrix()
        self.updateVoltages()
        self.updateEquationMatrix()

        return self.deltaMatrix, self.PEquationList, self.QEquationList

    #update the inverse jacobian matrix
    def updateInverseJacobian(self):
        numPQ = np.count_nonzero(self.busType == 'D')
        numPV = np.count_nonzero(self.busType == 'G')
        jacobian = np.zeros(((self.numBusses-1) + numPQ, (self.numBusses-1) + numPQ))
        jacobianSum = 0

        # H matrix
        # H goes from 2 to the number of busses
        for i in range(1, self.numBusses):
            for k in range(1, self.numBusses):
                if i != k:
                    jacobian[k - 1][i - 1] = (self.voltages[k] * self.voltages[i]) * (
                            self.admittanceReal[k][i] * np.sin(self.angles[k] - self.angles[i]) -
                            self.admittanceImag[k][i] * np.cos(self.angles[k] - self.angles[i])
                        )
                
                    jacobian[i - 1][i - 1] -= jacobian[k][i]

        # M matrix
        for i in range(numPV + 1, self.numBusses):
            for k in range(1, self.numBusses):
                if i != k:
                    jacobian[k - 1][i + numPQ - 1] = self.voltages[k] * (
                            self.admittanceReal[k][i] * np.cos(self.angles[k] - self.angles[i]) +
                            self.admittanceImag[k][i] * np.sin(self.angles[k] - self.angles[i])
                        )
                
                    jacobian[i + numPQ - 1][i + numPQ - 1] += (
                            jacobian[k - 1][i + numPQ - 1] +  2 * self.admittanceReal[k][k] * self.voltages[k]
                        )

        # N matrix
        for i in range(1, self.numBusses):
            for k in range(numPV + 1, self.numBusses):
                if i != k:
                    jacobian[k + numPQ - 1][i - 1] = (self.voltages[k] * self.voltages[i]) * (
                            np.negative(self.admittanceReal[k][i]) * np.cos(self.angles[k] - self.angles[i]) -
                            self.admittanceImag[k][i] * np.sin(self.angles[k] - self.angles[i])
                        )
                
                    jacobian[i - 1][i - 1] -= jacobian[k + numPQ - 1][i - 1]

        # L matrix
        for i in range(numPV + 1, self.numBusses):
            for k in range(numPV + 1, self.numBusses):
                if i != k:
                    jacobian[k + numPQ - 1][i + numPQ - 1] = self.voltages[k] * (
                            self.admittanceReal[k][i] * np.sin(self.angles[k] - self.angles[i]) -
                            self.admittanceImag[k][i] * np.cos(self.angles[k] - self.angles[i])
                        )

                    jacobianSum += jacobian[k + numPQ - 1][i + numPQ - 1] 
                    
            jacobian[i + numPQ - 1][i + numPQ - 1] = np.negative(2) * self.admittanceImag[k][k] * self.voltages[k] + jacobianSum

        #self.inverseJacobian = np.linalg.inv(jacobian) erroring out
        self.inverseJacobian = jacobian

    def updateDeltaMatrix(self):
        #check the order here, also should these be vertical matrices?
        self.deltaMatrix = np.dot( (-1 * self.inverseJacobian), self.implicitEquationList)

    def updateVoltages(self):
        #these won't be the same shape, I don't think. Will not be this simple
        #self.voltages += self.deltaMatrix
        #self.angles += self.deltaMatrix

        #for angle in deltaMatrix:
            #update angles

        #for voltage in deltaMatrix:
            #update voltages

        for i in range(1, self.numBusses):
            print(self.deltaMatrix[i])

    def updateEquationMatrix(self):
        
        self.PEquationList = np.zeros(self.numBusses)
        self.QEquationList = np.zeros(self.numBusses)

        self.implicitEquationList = []

        #create P and Q equation list
        for i in range(self.numBusses):
            for j in range(self.numBusses):

                PSum = (self.voltages[i] * self.voltages[j]) * (
                            self.admittanceReal[i][j]*np.cos(self.angles[i] - self.angles[j]) +
                            self.admittanceImag[i][j]*np.sin(self.angles[i] - self.angles[j])
                        )

                QSum = (self.voltages[i] * self.voltages[j]) * (
                            self.admittanceReal[i][j]*np.sin(self.angles[i] - self.angles[j]) -
                            self.admittanceImag[i][j]*np.cos(self.angles[i] - self.angles[j])
                        )
                    
                self.PEquationList[i] += PSum
                self.QEquationList[i] += QSum

            #add necessary P and Q equations to implicit equation list
                #SHOULD THIS BE A VERTICAL MATRIX?
            if(self.busType[i] == 'D'): #if PQ bus, both the p and q equation will be explicit
                self.implicitEquationList.append(self.PEquationList[i] - self.BusData['P MW'][i])
                self.implicitEquationList.append(self.QEquationList[i] - self.BusData['Q MVAr'][i])
            elif(self.busType[i] == 'G'): #if PV bus, only the p equation will be explicit
                self.implicitEquationList.append(self.PEquationList[i] - self.BusData['P MW'][i])

            #if slack bus, none of the equations are explicit
