import power_flow as pf

#test

#all that the constructor will do is initialize self variables, no processing
Solver = pf.PowerFlow(5, 0.1)

#this will udpate internal admittance matrix and create power flow equation matrix
Solver.readFromFile("notes_example.xlsx")

print(Solver.admittanceReal)
print(Solver.admittanceImag)


#while(True):
#    currentDelta, currentError = Solver.newtonRaphsonIteration()
#    
#    #checks for if NR should stop iterating, max iterations or solved
#    if(max(currentError) < tolerance):
#        print(currentDelta, currentError)
#        break
#    elif(Solver.currIterations > Solver.maxIterations):
#        print("Maximum iterations reached with no solutions")
#        break


#what do we actually need to output? Excel file, with what in it?
    #Network admittance matrix - already make that with Sam's function
    #Convergence history of bus results. Will need to store this
        #and update it every loop through NR
    #Final bus results. Will need to have a full equation matrix stored
        #as well as a equation matrix with only the implicit needed for NR
    # Line flows and MVA limit check... will need functions in class to do this as well

