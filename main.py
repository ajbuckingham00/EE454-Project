#This file controls the operation of a power flow solving program
#The functions of the program are contained in a class in a separate file called power_flow.py
#This main.py file will setup the class with input data, then manage the iterations of the newton raphson technique
#   which is used to solve the power flow

import power_flow as pf
import numpy as np

#create a new power flow solver, with arguments for maximum number of iterations and mismatch tolerance before stopping iterations
Solver = pf.PowerFlow(5, 0.001)

#update power flow solver with input data from excel file given in argument
Solver.readFromFile("system_basecase.xlsx")

i = 1 #index counter for tracking number of iterations

#loop will perform newton raphson iterations until a solution converges or the maximum number of iterations is reached
while(True): 
    
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Iteration ", i)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    i += 1
    
    currentMismatches = Solver.newtonRaphsonIteration() #get the current mismatch values to see if they are below the tolerance

    #if the current mismatches are below the tolerance, stop iterating because a solution is reached
    #otherwise, keep iterating until maximum iterations are reached or a solution is found
    if(max(np.absolute(currentMismatches)) < Solver.tolerance):
        
        print("Success!")
        print()
        break
    
    elif(Solver.currIterations > Solver.maxIterations):
        print("Maximum iterations reached with no solutions")
        break


#store output data into excel file with name given in argument. 
Solver.output("output_base.xlsx")
