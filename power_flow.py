import numpy as np
import pandas as pd

class PowerFlow():

    def __init__(self, maxIterations, tolerance):

        self.maxIterations = maxIterations
        self.currIterations = 0
        
        self.tolerance = tolerance

        self.admittanceReal = 0
        self.admittanceImag = 0
        
        self.equationMatrix = 0
        self.implicitEquationMatrix = 0
        
        self.deltaMatrix = 0
        self.inverseJacobian = 0

        self.voltages = 0
        self.angles = 0
        self.busType = 0

        self.convergenceHistory = 0

    #pulls the 
    def readFromFile(self, fileName):
        BusData = pd.read_excel(fileName, sheet_name = 0)
        LineData = pd.read_excel(fileName, sheet_name = 1)

        BusNumber = len(BusData) # counts the number of buses
        LineNumber = len(LineData) # counts the number of connections

        print(LineData['From'].loc[LineData.index[0]])

        # creates a zero matrix based on the number of buses

        self.admittanceReal = np.zeros((BusNumber, BusNumber)) 
        self.admittanceImag = np.zeros((BusNumber, BusNumber))

        # fills in the line connection values to the matrix
    
        for i in range(LineNumber):
            ind1 = LineData['From'].loc[LineData.index[i]]
            ind2 = LineData['To'].loc[LineData.index[i]]
            real = LineData['Rtotal, p.u.'].loc[LineData.index[i]]
            imag = LineData['Xtotal, p.u.'].loc[LineData.index[i]]

            if imag == 0:
                new_imag = 0
            else:
                new_imag = np.negative(1/imag)
            if real == 0:
                new_real = 0
            else:
                new_real = 1/real


            self.admittanceImag[ind1 - 1, ind2 - 1] = np.negative(new_imag)
            self.admittanceImag[ind2 - 1, ind1 - 1] = np.negative(new_imag)
            self.admittanceReal[ind1 - 1, ind2 - 1] = np.negative(new_real)
            self.admittanceReal[ind2 - 1, ind1 - 1] = np.negative(new_real)
        
            # Starts adding in the diagonal numbers
            
            for i in range(BusNumber):
                print("i: ", i)
                Real_rowtotal = 0
                Imag_rowtotal = 0
                b = LineData['Btotal, p.u.'].loc[LineData.index[i]]
                
                if b == 0:
                    new_b = 0
                else:
                    new_b = b / 2
                    
                
                for j in range(BusNumber):
                    Real_rowtotal += np.negative(self.admittanceReal[i][j])
                    Imag_rowtotal += np.negative(self.admittanceImag[i][j])
                    
                self.admittanceReal[i, i] = Real_rowtotal
                self.admittanceImag[i, i] = Imag_rowtotal + new_b

        print(self.admittanceReal, self.admittanceImag)

        #add functionality here... pull voltages from busdata into voltages matrix for NR
        self.voltages = [[]]
        self.angles = [[]]
        self.busType = [[]]

        #update the output values of the functions to be solved with the initial flat case as input
        self.updateEquationMatrix()

    def newtonRaphsonIteration(self):
        #update all these values
        self.currIterations += 1
        
        self.updateinverseJacobian()
        self.updateDeltaMatrix()
        self.updateVoltages()
        self.updateEquationMatrix()

        return self.deltaMatrix, self.equationMatrix

    def updateInverseJacobian(self):
        return
        #for loop through each term and build up Jacobian according to power flow slide 100-104
        #jacobian = ...
        #self.inverseJacobian = np.linalg.inv(jacobian)

    def updateDeltaMatrix(self):
        return
        #self.deltaMatrix = np.dot( (-1 * self.inverseJacobian), self.equationMatrix)

    def updateVoltages(self):
        return
        #these won't be the same shape, I don't think. Will not be this simple
        #self.voltages += self.deltaMatrix

    def updateEquationMatrix(self):
        return
        #for loop through matrix and update each equation according to power flow slides
        #if the equation has a known P or Q, then also add it to the implicit equations list
        #will have to look at the type of bus to do this, I think?
        #also seems like we will end up with too many implicit equations, how do we pick which ones to use?

        #self.implicitEquationMatrix = ...

        #for i in range():
        #    for j in range():
        #        val = #value for that cell
        #        if(busType[][] == ...):
                    #if the bus type indicates this will be an explicit equation, also add it to explicit list

        
        #will have to figure out which ones are implicit and which are explicit, only put implicit equations into this matrix
        #how will we do that? look at the busdata and see where we know P and Q, cross reference with list?
        
        #self.equationMatrix = ... #some subset of the full equation matrix with only the implicit equations needed for NR
    
        
