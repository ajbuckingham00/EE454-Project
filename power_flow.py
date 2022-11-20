class PowerFlow():

    def __init__(maxIterations, tolerance):

        self.maxIterations = maxIternations
        self.currIterations = 0
        
        self.tolerance = tolerance

        self.admittanceReal
        self.admittanceImag
        
        self.equationMatrix
        self.deltaMatrix
        self.voltages
        self.inverseJacobian

    #pulls the 
    def readFromFile(fileName):
        BusData = pd.read_excel(fileName, sheet_name = 0)
        LineData = pd.read_excel(fileName, sheet_name = 1)

        BusNumber = len(BusData) # counts the number of buses
        LineNumber = len(LineData) # counts the number of connections

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

        #add functionality here... pull voltages from busdata into voltages matrix for NR
        self.voltages = ...

        #update the output values of the functions to be solved with the initial flat case as input
        self.updateEquationMatrix()

    def newtonRaphsonIteration:
        #update all these values
        self.currIterations += 1
        
        self.inverseJacobian = self.updateinverseJacobian()
        self.deltaMatrix = self.updateDeltaMatrix()
        self.voltages = self.updateVoltages()
        self.equationMatrix = self.updateEquationMatrix()

        return self.deltaMatrix, self.equationMatrix

    def updateInverseJacobian():
        #for loop through each term and build up Jacobian according to power flow slide 100-104
        jacobian = ...
        return np.linalg.inv(jacobian)

    def updateDeltaMatrix():
        return np.dot( (-1 * self.inverseJacobian), self.equationMatrix)

    def updateVoltages():
        #these won't be the same shape, I don't think. Will not be this simple
        self.voltages += self.deltaMatrix

    def updateEquationMatrix():
        #for loop through matrix and update each equation according to power flow slides
        fullEquationMatrix = ...

        
        #will have to figure out which ones are implicit and which are explicit, only put implicit equations into this matrix
        #how will we do that? look at the busdata and see where we know P and Q, cross reference with list?
        
        self.equationMatrix = ... #some subset of the full equation matrix with only the implicit equations
    
        
