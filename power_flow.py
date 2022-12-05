import numpy as np
import pandas as pd
from openpyxl import Workbook
import math


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
        self.numPQ = 0
        self.numPV = 0

        self.busMap = 0

        self.maxMismatchP = []
        self.maxMismatchPLocations = []
        self.maxMismatchQ = []
        self.maxMismatchQLocations = []

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
            node1Index = np.where(self.busMap == (node1 - 1))[0][0]
            node2Index = np.where(self.busMap == (node2 - 1))[0][0]

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

        print("Voltages: ")
        print(self.voltages)
        print()

        print("Angles: ")
        print(self.angles)
        print()

    def getBusType(self):
        self.busType = np.array(self.BusData['Type'])
        self.numPQ = np.count_nonzero(self.busType == 'D')
        self.numPV = np.count_nonzero(self.busType == 'G')

    def newtonRaphsonIteration(self):
        #update all these values
        self.currIterations += 1
        
        self.updateInverseJacobian()
        self.updateDeltaList()
        self.updateVoltages()
        self.updateEquationLists()

        return self.deltaList, self.implicitEquationList

    #update the inverse jacobian matrix
    def updateInverseJacobian(self):
        
        jacobian = np.zeros( ( (self.numBusses-1) + self.numPQ, (self.numBusses-1) + self.numPQ), dtype = float)

        jacobian = self.buildHMatrix(jacobian)
        jacobian = self.buildMMatrix(jacobian)
        jacobian = self.buildNMatrix(jacobian)
        jacobian = self.buildLMatrix(jacobian)

        # print("Jacobian:")
        # print(jacobian)
        print()
        self.inverseJacobian = np.linalg.inv(jacobian)
        
    def buildHMatrix(self, jacobian):

        jacobianRow = 0 #the starting row of the jacobian for region H
        
        for k in range(2, self.numBusses + 1): #k goes from 2 to N
            k -= 1 #to use k for matrix indices that start at 0
            jacobianCol = 0 #the starting column of the jacobian for region H
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
        
        for k in range(2, self.numBusses + 1): #k goes from 2 to N
            k -= 1 #to use k for matrix indices that start at 0
            jacobianCol = self.numBusses - 1 #the starting column of the jacobian for region M
            for i in range((self.numPV + 1) + 1, self.numBusses + 1): #i goes from m+1 to N, m is number of PV buses + 1
                i -= 1 #to use i for matrix indices that start at 0

                if i != k:
                    jacobian[jacobianRow][jacobianCol] = self.voltages[k] * (
                            self.admittanceReal[k][i] * np.cos(self.angles[k] - self.angles[i]) +
                            self.admittanceImag[k][i] * np.sin(self.angles[k] - self.angles[i])
                        )
                    
                    
                else:
                    jacobian[jacobianRow][jacobianCol] += 2 * self.admittanceReal[k][k] * self.voltages[k]
                    for l in range(1, self.numBusses + 1):
                        l -= 1 #to use l for matrix indices that start at 0
                        if l != k:
                            jacobian[jacobianRow][jacobianCol] += self.voltages[l] * (
                                    self.admittanceReal[k][l] * np.cos(self.angles[k] - self.angles[l]) +
                                    self.admittanceImag[k][l] * np.sin(self.angles[k] - self.angles[l])
                                )
                jacobianCol += 1
            jacobianRow += 1

        return jacobian

    def buildNMatrix(self, jacobian):

        jacobianRow = self.numBusses - 1 #the starting row of the jacobian for region N
        
        for k in range((self.numPV + 1) + 1, self.numBusses + 1): #k goes from m+1 to N, m is number of PV buses + 1
            k -= 1 #to use k for matrix indices that start at 0
            jacobianCol = 0 #the starting column of the jacobian for region N
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
        
        for k in range((self.numPV + 1) + 1, self.numBusses + 1): #k goes from m+1 to N, m is number of PV buses + 1
            k -= 1 #to use k for matrix indices that start at 0
            jacobianCol = self.numBusses - 1 #the starting column of the jacobian for region L
            for i in range((self.numPV + 1) + 1, self.numBusses + 1): #k goes from m+1 to N, m is number of PV buses + 1
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
        print(self.implicitEquationList)
        print("vstack")
        print( np.vstack(self.implicitEquationList))
        print()
        
        self.deltaList = np.dot( (-1 * self.inverseJacobian), self.implicitEquationList)
        print("Delta List:")
        print(self.deltaList)
        print()

    def updateVoltages(self):
        #delta list will have theta2 to thetaN, then vm+1 to vN

        angleIndex = 1 #index of first angle to be modified by delta list
        for i in range(self.numBusses - 1): #pull angles out of delta list and add them to angles list
            self.angles[angleIndex] += self.deltaList[i]
            angleIndex += 1

        voltageIndex = self.numPV + 1 #index of first voltage to be modified by delta list
        for i in range(self.numBusses - 1, len(self.deltaList)): #pull voltages out of delta list and add them to voltages list
            self.voltages[voltageIndex] += self.deltaList[i] 
            voltageIndex += 1

        print("Voltages:")
        print(self.voltages)
        print()

        print("Angles:")
        print(self.angles)
        print()

        print(self.busMap)
        print()

    def updateEquationLists(self):
        
        self.PEquationList = np.zeros(self.numBusses)
        self.QEquationList = np.zeros(self.numBusses)

        self.implicitEquationList = [] 
        #create P and Q equation list
        for i in range(self.numBusses):
            for j in range(self.numBusses):
                PSum = (self.voltages[i] * self.voltages[j]) * (
                            self.admittanceReal[i][j] * np.cos(self.angles[i] - self.angles[j]) +
                            self.admittanceImag[i][j] * np.sin(self.angles[i] - self.angles[j])
                        )

                QSum = (self.voltages[i] * self.voltages[j]) * (
                            self.admittanceReal[i][j] * np.sin(self.angles[i] - self.angles[j]) -
                            self.admittanceImag[i][j] * np.cos(self.angles[i] - self.angles[j])
                        )
                
                self.PEquationList[i] += PSum
                self.QEquationList[i] += QSum
            print("power: ", self.PEquationList[i])
            #print("reactive: ", self.QEquationList[i])
            print()

        maximumPMismatch = 0
        currPMismatch = 0
        maximumPLocation = 0
            
        for i in range(len(self.PEquationList)): #get each item in the P list, see if it is implicit or not and add it if so
            if(self.busType[i] != 'S'): #both PQ and PV buses have an explicit P equation
                currPMismatch = self.PEquationList[i] - ((self.BusData['P Gen'][self.busMap[i]] / 100) - (self.BusData['P MW'][self.busMap[i]] / 100))
                if currPMismatch > maximumPMismatch:
                    maximumPMismatch = currPMismatch
                    maximumPLocation = i
                    
                self.implicitEquationList.append([currPMismatch])

        maximumQMismatch = 0
        currQMismatch = 0
        maximumQLocation = 0

        for i in range(len(self.QEquationList)): #get each item in the Q list, see if it is implicit or not and add it if so
            if(self.busType[i] == 'D'): #only PQ buses have an explicit Q equation
                currQMismatch = self.QEquationList[i] + (self.BusData['Q MVAr'][self.busMap[i]] / 100)
                if currQMismatch > maximumQMismatch:
                    maximumQMismatch = currQMismatch
                    maximumQLocation = i
                    
                self.implicitEquationList.append([currQMismatch])

        self.maxMismatchP.append(maximumPMismatch)
        self.maxMismatchQ.append(maximumQMismatch)
        self.maxMismatchPLocations.append(maximumPLocation)
        self.maxMismatchQLocations.append(maximumQLocation)

        #if slack bus, none of the equations are explicit
        print("implicit equations...")
        print(np.vstack(self.implicitEquationList))

    def output(self, filename, mismatch):
        wb = Workbook()
        sheet1 = wb.create_sheet("List of Voltages", 0)
        sheet2 = wb.create_sheet("Power Produce at Gen.", 1)
        sheet3 = wb.create_sheet("Power Flow at Line", 2)
        sheet4 = wb.create_sheet("Line Check", 3)
        sheet5 = wb.create_sheet("Voltage Check", 4)
        sheet6 = wb.create_sheet("Convergence Record", 5)

        #build sheet1
        for i in range(self.numBusses):
            voltageStatement = "Bus ", i + 1, "Voltage: ", self.voltages[np.where(self.busMap == i)[0][0]], "Angle: ", self.angles[np.where(self.busMap == i)[0][0]] * (180 / math.pi)
            sheet1.append(voltageStatement)

        #build sheet2
        for i in range(len(self.PEquationList)):
            flowStatement = "Bus ", i + 1, "P MW: ", self.PEquationList[i] * 100, "Q MVAr: ", self.QEquationList[i] * 100 
            sheet2.append(flowStatement)

        #build sheet3 and sheet4
        for i in range(self.numConnections):
            node1 = self.LineData['From'][i]
            node2 = self.LineData['To'][i]
            
            node1Index = np.where(self.busMap == node1 - 1)[0][0]
            node2Index = np.where(self.busMap == node2 - 1)[0][0]

            print(node1Index)
            print(node2Index)

            realPower = (self.voltages[node1Index] * self.voltages[node2Index]) * (
                            self.admittanceReal[node1Index][node2Index] * np.cos(self.angles[node1Index] - self.angles[node2Index]) + 
                            self.admittanceImag[node1Index][node2Index] * np.sin(self.angles[node1Index] - self.angles[node2Index])
                        )
            
            reactivePower = (self.voltages[node1Index] * self.voltages[np.where(self.busMap == node2 - 1)[0][0]]) * (
                                self.admittanceReal[node1Index][node2Index] * np.sin(self.angles[node1Index] - self.angles[node2Index]) -
                                self.admittanceImag[node1Index][node2Index] * np.cos(self.angles[node1Index] - self.angles[node2Index])
                            )
            MVA = math.sqrt((100 * realPower) ** 2 + (100 * reactivePower) ** 2)

            lineStatement = "Line:", node1, "to", node2, "P MW:", realPower * 100, "Q MVAr:", reactivePower * 100, "S MVA:", MVA
            sheet3.append(lineStatement)

            if  MVA > self.LineData['Fmax, MVA'][i]:
                lineData = "Line ", node1, " to Line ", node2, " has exceeded normal operation limit"
                sheet4.append(lineData)
                
            else:
                lineData = "Line ", node1, " to Line ", node2, " has not exceeded normal operation limit"
                sheet4.append(lineData)

        #build sheet 5
        for i in range(self.numBusses):
            if self.voltages[np.where(self.busMap == i)[0][0]] > 1.05 or self.voltages[np.where(self.busMap == i)[0][0]] < 0.95:
                voltageData = "Bus ", i + 1, "has exceeded normal operating limits"
                sheet5.append(voltageData)
            else:
                voltageData = "Bus ", i + 1, "has not exceeded normal operating limits"
                sheet5.append(voltageData)

        #build sheet6
        sheet6.append(["Iteration", "Maximum P Mismatch", "Location of P Mismatch", "Maximum Q Mismatch", "Location of Q Mismatch"])

        #trim extra mismatch data from mismatches before first newton raphson iteration
        self.maxMismatchP = self.maxMismatchP[1:len(self.maxMismatchP)]
        self.maxMismatchPLocations = self.maxMismatchPLocations[1:len(self.maxMismatchPLocations)]
        self.maxMismatchQ = self.maxMismatchQ[1:len(self.maxMismatchQ)]
        self.maxMismatchQLocations = self.maxMismatchQLocations[1:len(self.maxMismatchQLocations)]

        for i in range(len(self.maxMismatchP)):
            sheet6.append([i + 1, self.maxMismatchP[i], self.maxMismatchPLocations[i], self.maxMismatchQ[i], self.maxMismatchQLocations[i]])

        wb.save(filename)
