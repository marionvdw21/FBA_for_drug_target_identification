#from GetEssRxn2 import *
#from convert2NALM import *
import cobra.io as co
import random
#from GeneInteraction3 import *
#from findTarget4 import * 
import cobra
from convert2MM import convertMM2
from convert2NALM6 import *
#from GOT1 import etc

def makeLayover(): 
    #~~~CHOICE OF MODEL~~~~
    #model = convertToNalm()[0]
    model = convertToNalm6('expression_Recon_ordered2.txt')[0]
    #model = convertMM2(convertToNalm6('expression_Recon_deordered2.txt')[0])
    model = co.read_sbml_model('Recon3.xml')
    #model.reactions.EX_o2s_e.upper_bound = 20
    #model.reactions.EX_o2s_e.lower_bound = -20
    #model.reactions.EX_o2_e.upper_bound = 0
    #model.reactions.EX_o2_e.lower_bound = 0

    color_list = ['800000', '9A6324', '808000', '469990', '000075', '000000', 'e6194B', 'f58231', '5DFDCB', 'ffe119', '3cb44b', '42d4f4', '4363d8', '911eb4', 'f032e6', 'a9a9a9', 'fabed4', 'ffd8b1', 'aaffc3', 'dcbeff'];  
    #color = color_list[random.randint(0, 19)]
    color = 'A5CC6B'
    #maybe fo an else statement for errors. 
    rxnList = [] # List I use t store the lines i'm going to write in 'layover.txt'
    #rxnListCheck = [] # List to check that I don't write the same rxn 2 times
    header = ''
    header += 'name' + '\t' + 'reactionIdentifier' + '\t' + 'lineWidth' + '\t' + 'color'
    rxnList.append(header)
    #geneList = faceCBD()

    #~~~CHOICE OF FLUX~~~
    #myFlux = findTarget4()[1][1]
    myFlux = model.optimize().fluxes
    #myFlux = cobra.flux_analysis.pfba(model=model).fluxes
    #WTsolution = model1.optimize()
    #myFlux = cobra.flux_analysis.moma(model, WTsolution).fluxes
    myRxn = myFlux.index
    
    for i in range(len(myFlux)): 
        if myFlux[i] != 0.0:
            myString= ''
            reactionName = ''
            reactionName += 'R_' + str(myRxn[i])
            LW = abs(float(myFlux[i])) * 0.005
            myString += '\t' + reactionName + '\t' + str(5) + '\t' + color #I could adjust thickness of line depending on flux
            rxnList.append(myString)
   
    OutputFileName = 'layover_%(someGene)s.txt' %{'someGene' : 'nalm6model'} 
    with open(OutputFileName, 'w') as file: 
        for k in range(len(rxnList)): 
            file.write(rxnList[k])
            file.write('\n')

    return 
makeLayover()
