def makeLayover(myModel): 
    import random

    model = myModel
   

    color_list = ['800000', '9A6324', '808000', '469990', '000075', '000000', 'e6194B', 'f58231', '5DFDCB', 'ffe119', '3cb44b', '42d4f4', '4363d8', '911eb4', 'f032e6', 'a9a9a9', 'fabed4', 'ffd8b1', 'aaffc3', 'dcbeff'];  
    color = color_list[random.randint(0, 19)]

    rxnList = [] # List I use t store the lines i'm going to write in 'layover.txt'
    header = ''
    header += 'name' + '\t' + 'reactionIdentifier' + '\t' + 'lineWidth' + '\t' + 'color'
    rxnList.append(header)


    #~~~CHOICE OF FLUX~~~
    myFlux = model.optimize().fluxes
    myRxn = myFlux.index
    
    for i in range(len(myFlux)): 
        if myFlux[i] != 0.0:
            myString= ''
            reactionName = ''
            reactionName += 'R_' + str(myRxn[i])
            LW = abs(float(myFlux[i])) * 0.005
            myString += '\t' + reactionName + '\t' + str(5) + '\t' + color #I could adjust thickness of line depending on flux
            rxnList.append(myString)
   
    OutputFileName = 'layover_%(someGene)s.txt' %{'someGene' : 'myModel'} 
    with open(OutputFileName, 'w') as file: 
        for k in range(len(rxnList)): 
            file.write(rxnList[k])
            file.write('\n')

    return 

if __name__ == '__main__': 
    import convertMyMedia
    import convertMyModel
    import cobra.io as co
    myModel = convertMyMedia.convertMyMedia(convertMyModel.convertMyModel('expression_Recon_ordered2.txt', co.read_sbml_model('Recon3.xml'))[0], 'minimal_media2.txt')
    makeLayover(myModel)
