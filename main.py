import power_flow as pf

#all that the constructor will do is initialize self variables, no processing
Solver = new PowerFlow(maxIterations, tolerance)

#this will udpate internal admittance matrix and create power flow equation matrix
Solver.readFromFile()


while(True):
    currentDelta, currentError = Solver.newtonRaphsonIteration()
    
    #checks for if NR should stop iterating, max iterations or solved
    if(Solver.currIterations > Solver.maxIterations):
        #raise error, or something. did not converge
        break
    if(min(currentError) < tolerance):
        #update output. solution is found
        break


#what do we actually need to output? Excel file, with what in it?
    #Network admittance matrix - already make that with Sam's function
    #Convergence history of bus results. Will need to store this
        #and update it every loop through NR
    #Final bus results. Will need to have a full equation matrix stored
        #as well as a equation matrix with only the implicit needed for NR
    # Line flows and MVA limit check... will need functions in class to do this as well

