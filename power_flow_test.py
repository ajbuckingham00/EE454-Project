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
        
        self.deltaLista = 0
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
        self.updateEquationLists()

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

        print("Real Admittance:")
        print(self.admittanceReal)
        print()
        print("Imaginary Admittance:")
        print(self.admittanceImag)
        print()
        
    def getVoltages(self): #voltages AND angles
        self.voltages = np.array(self.BusData['V Set'], dtype = float)
        self.angles = np.zeros_like(self.voltages, dtype = float)

    def getBusType(self):
        self.busType = np.array(self.BusData['Type'])

    def newtonRaphsonIteration(self):
        #update all these values
        self.currIterations += 1
        
        self.updateInverseJacobian()
        self.updateDeltaList()
        self.updateVoltages()
        self.updateEquationLists()

        return self.deltaList, self.implicitEquationList

    #update the inverse jacobian matrix

    #SCRAP AND START OVER... Hik can't mean the actual position within H, has to just be i and k in the equation and its linked to position by the equality with the derivative
    #this is different technique than both tried so far
    def updateInverseJacobian(self):
        numPQ = np.count_nonzero(self.busType == 'D')
        numPV = np.count_nonzero(self.busType == 'G')
        jacobian = np.zeros(((self.numBusses-1) + numPQ, (self.numBusses-1) + numPQ), dtype = float)

        jacobian = self.buildHMatrix(jacobian)
        jacobian = self.buildMMatrix(jacobian)
        jacobian = self.buildNMatrix(jacobian)
        jacobian = self.buildLMatrix(jacobian)

        print("Jacobian:")
        print(jacobian)
        print()
        self.inverseJacobian = np.linalg.inv(jacobian)
        
    def buildHMatrix(self, jacobian):

        jacobianRow = 0 #the starting row of the jacobian for region H
        jacobianCol = 0 #the starting column of the jacobian for region H
        
        for k in range(2, self.numBusses + 1): #k goes from 2 to N
            k -= 1 #to use k for matrix indices that start at 0
            for i in range(2, self.numBusses + 1): #i goes from 2 to N
                i -= 1 #to use i for matrix indices that start at 0

                if i != k:
                    jacobian[jacobianRow][jacobianCol] = self.voltages[k]*self.voltages[i] * (
                            self.admittanceReal[k][i] * np.sin(self.angles[k] - self.angles[i]) -
                            self.admittanceImag[k][i] * np.cos(self.angles[k] - self.angles[i])
                        )
                    
                else:
                    for l in range(1, self.numBusses + 1):
                        l -= 1 #to use l for matrix indices that start at 0
                        if l != k:
                            jacobian[jacobianRow][jacobianCol] += self.voltages[k]*self.voltages[l] * (
                                    -1.0 * self.admittanceReal[k][l] * np.sin(self.angles[k] - self.angles[l]) +
                                    self.admittanceImag[k][l] * np.cos(self.angles[k] - self.angles[l])
                                )
                jacobianCol += 1
            jacobianRow += 1
            
        return jacobian

    def buildMMatrix(self, jacobian):

        jacobianRow = 0 #the starting row of the jacobian for region M
        jacobianCol = self.numBusses - 1 #the starting column of the jacobian for region M
        
        for k in range(2, self.numBusses + 1): #k goes from 2 to N
            k -= 1 #to use k for matrix indices that start at 0
            for i in range((np.count_nonzero(self.busType == 'G') + 1) + 1, self.numBusses + 1): #i goes from m+1 to N, m is number of PV buses + 1
                i -= 1 #to use i for matrix indices that start at 0

                if i != k:
                    jacobian[jacobianRow][jacobianCol] = self.voltages[k] * (
                            self.admittanceReal[k][i] * np.cos(self.angles[k] - self.angles[i]) +
                            self.admittanceImag[k][i] * np.sin(self.angles[k] - self.angles[i])
                        )
                    
                else:
                    for l in range(1, self.numBusses + 1):
                        l -= 1 #to use l for matrix indices that start at 0
                        if l != k:
                            jacobian[jacobianRow][jacobianCol] += self.voltages[l] * (
                                    self.admittanceReal[k][l] * np.cos(self.angles[k] - self.angles[l]) +
                                    self.admittanceImag[k][l] * np.sin(self.angles[k] - self.angles[l])
                                ) + 2 * self.admittanceReal[k][k] * self.voltages[k]
                jacobianCol += 1
            jacobianRow += 1

        return jacobian

    def buildNMatrix(self, jacobian):

        jacobianRow = self.numBusses - 1 #the starting row of the jacobian for region N
        jacobianCol = 0 #the starting column of the jacobian for region N
        
        for k in range((np.count_nonzero(self.busType == 'G') + 1) + 1, self.numBusses + 1): #k goes from m+1 to N, m is number of PV buses + 1
            k -= 1 #to use k for matrix indices that start at 0
            for i in range(2, self.numBusses + 1): #i goes 2 to N
                i -= 1 #to use i for matrix indices that start at 0

                if i != k:
                    jacobian[jacobianRow][jacobianCol] = self.voltages[k]*self.voltages[i] * (
                            -1.0 * self.admittanceReal[k][i] * np.cos(self.angles[k] - self.angles[i]) +
                            -1.0 * self.admittanceImag[k][i] * np.sin(self.angles[k] - self.angles[i])
                        )
                    
                else:
                    for l in range(1, self.numBusses + 1):
                        l -= 1 #to use l for matrix indices that start at 0
                        if l != k:
                            jacobian[jacobianRow][jacobianCol] += self.voltages[k]*self.voltages[l] * (
                                    self.admittanceReal[k][l] * np.cos(self.angles[k] - self.angles[l]) +
                                    self.admittanceImag[k][l] * np.sin(self.angles[k] - self.angles[l])
                                )
                jacobianCol += 1
            jacobianRow += 1

        return jacobian

    def buildLMatrix(self, jacobian):

        jacobianRow = self.numBusses - 1 #the starting row of the jacobian for region L
        jacobianCol = self.numBusses - 1 #the starting column of the jacobian for region L
        
        for k in range((np.count_nonzero(self.busType == 'G') + 1) + 1, self.numBusses + 1): #k goes from m+1 to N, m is number of PV buses + 1
            k -= 1 #to use k for matrix indices that start at 0
            for i in range((np.count_nonzero(self.busType == 'G') + 1) + 1, self.numBusses + 1): #k goes from m+1 to N, m is number of PV buses + 1
                i -= 1 #to use i for matrix indices that start at 0

                if i != k:
                    jacobian[jacobianRow][jacobianCol] = self.voltages[k] * (
                            self.admittanceReal[k][i] * np.sin(self.angles[k] - self.angles[i]) +
                            -1.0 * self.admittanceImag[k][i] * np.cos(self.angles[k] - self.angles[i])
                        )
                    
                else:
                    jacobian[jacobianRow][jacobianCol] += -2 * self.admittanceImag[k][k]*self.voltages[k]
                    for l in range(1, self.numBusses + 1):
                        l -= 1 #to use l for matrix indices that start at 0
                        if l != k:
                            jacobian[jacobianRow][jacobianCol] += self.voltages[l] * (
                                    self.admittanceReal[k][l] * np.sin(self.angles[k] - self.angles[l]) -
                                    self.admittanceImag[k][l] * np.cos(self.angles[k] - self.angles[l])
                                )
                jacobianCol += 1
            jacobianRow += 1

        return jacobian
                    
                

    def updateDeltaList(self):
        #check the order here, also should these be vertical matrices?
        self.deltaList = np.dot( (-1 * self.inverseJacobian), self.implicitEquationList)
        print("Delta List:")
        print(self.deltaList)
        print()

    def updateVoltages(self):
        #delta list will have theta2 to thetaN, then vm+1 to vN

        angleIndex = 1 #index of first angle to be modified by delta list
        for i in range(1, self.numBusses): #pull angles out of delta list and add them to angles list
            self.angles[angleIndex] += self.deltaList[i][0]
            angleIndex += 1

        voltageIndex = np.count_nonzero(self.busType == 'G') + 1 #index of first voltage to be modified by delta list
        for i in range(self.numBusses - 1, self.numBusses + np.count_nonzero(self.busType == 'G')): #pull voltages out of delta list and add them to voltages list
            self.voltages[voltageIndex] += self.deltaList[i][0]
            voltageIndex += 1

        print("Voltages:")
        print(self.voltages)
        print()

        print("Angles:")
        print(self.angles)
        print()

    def updateEquationLists(self):
        
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

        for i in range(0, len(self.PEquationList)): #get each item in the P list, see if it is implicit or not and add it if so
            if(self.busType[i] != 'S'): #both PQ and PV buses have an explicit P equation
                self.implicitEquationList.append([self.PEquationList[i] - self.BusData['P MW'][i]])

        for i in range(0, len(self.QEquationList)): #get each item in the Q list, see if it is implicit or not and add it if so
            if(self.busType[i] == 'D'): #only PQ buses have an explicit Q equation
                self.implicitEquationList.append([self.QEquationList[i] - self.BusData['Q MVAr'][i]])

        #if slack bus, none of the equations are explicit
        print("implicit equations...")
        print(self.implicitEquationList)
