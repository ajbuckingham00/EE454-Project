#!/usr/bin/env python
# coding: utf-8

# In[5]:


import pandas as pd
import numpy as np

BusData = pd.read_excel("system_basecase.xlsx", sheet_name = 0)
LineData = pd.read_excel("system_basecase.xlsx", sheet_name = 1)

BusNumber = len(BusData) # counts the number of buses
LineNumber = len(LineData) # counts the number of connections

# creates a zero matrix based on the number of buses

RealMatrix = np.zeros((BusNumber, BusNumber)) 
ImagMatrix = np.zeros((BusNumber, BusNumber))




# In[3]:


# This function will take in Line Data and the number of rows of data 
# and create a matrix of real and imaginary admittance matrix.
# The function will look at which bus is connected and place the admittance into the matrix
# The function will also add all that connection and the shunt admittance and place them in
# the proper diagonal spot. 

def fill_Matrix(LineData, LineNumber):
    
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


        ImagMatrix[ind1 - 1, ind2 - 1] = np.negative(new_imag)
        ImagMatrix[ind2 - 1, ind1 - 1] = np.negative(new_imag)
        RealMatrix[ind1 - 1, ind2 - 1] = np.negative(new_real)
        RealMatrix[ind2 - 1, ind1 - 1] = np.negative(new_real)
    
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
            Real_rowtotal += np.negative(RealMatrix[i][j])
            Imag_rowtotal += np.negative(ImagMatrix[i][j])
            
        RealMatrix[i, i] = Real_rowtotal
        ImagMatrix[i, i] = Imag_rowtotal + new_b
        
    return RealMatrix, ImagMatrix # returns the finish matrix for both real and imaginary

#pull voltages from excel file, return matrix of voltages at each bus and 
def getVoltages(BusData, BusNumber) {
    #...
    #return k length list matrix of voltages, k length list of angles
}

#
def createPowerEquations(realMatrix, ImagMatrix, voltageMatrix, angleMatrix) {
    #return k length list of power flow equations, one for P, one for Q
    #
    #for.. k
    #   powerFlowEquationList[k] = 0
    #   for... i
    #   powerFlowEquationList[k] += voltageMatrix[k]voltageMatrix[i]*......
}

#
def newtonRaphson(powerFlowEquationList) {
    #for each equation in powerFlowEquationList...
    #   pBase = plug in base value for changing variable and solve equation
    #   pDiff = plug in base value + delta value for changing variable and solve equation
    #   derivative = (pDiff-pBase)/delta value

    #store derivatives in jacobian matrix

    #invert jacobian matrix

    #deltatheta, deltaV = multiply inverse jacobian by -(powerFlowEquations)
    
    #return powerFlowEquationList + deltatheta, deltav, list of difference in theta and v values
}


# calls the function to fill in the matrix with data

RealMatrix, ImagMatrix = fill_Matrix(LineData, LineNumber)

#call all functions up until NR once

while(diffV, difftheta > ...) {
    powerFlowEquationList, diffV, difftheta = newtonRaphson(powerFlowEquationList)
}
