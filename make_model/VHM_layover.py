def makeLayover(myModel): 
    '''
    U S A G E : 
    This function is used to create a text file compatible with the layover option of the Virtual Metabolic Human website (https://www.vmh.life/#reconmap). 
    This option allows the user to visualize which reactions in the metabolic model carry flux and which do not. 
    By creating an account >  logging in > 'overlays' > 'add overlay' > the user will be prompted to load a text file. 
    This makeLayover() function makes this text file which represents the flux passing through each reaction when a specific model 
    (inputed as myModel) is optimized. 

    I N P U T : 
    --> myModel : the model of which the layover needs to be created. 

    O U T P U T : 
    --> myModel_layover.txt : a text file containing the information for the layover (reaction ID, lineweight and color). 

    '''


    import random
    model = myModel
   

    color_list = ['800000', '9A6324', '808000', '469990', '000075', '000000', 'e6194B', 'f58231', '5DFDCB', 'ffe119', '3cb44b', '42d4f4', '4363d8', '911eb4', 'f032e6', 'a9a9a9', 'fabed4', 'ffd8b1', 'aaffc3', 'dcbeff'];  
    color = color_list[random.randint(0, 19)]

    rxnList = [] 
    header = ''
    header += 'name' + '\t' + 'reactionIdentifier' + '\t' + 'lineWidth' + '\t' + 'color'
    rxnList.append(header)


    myFlux = model.optimize().fluxes
    myRxn = myFlux.index
    
    for i in range(len(myFlux)): 
        if myFlux[i] != 0.0:
            myString= ''
            reactionName = ''
            reactionName += 'R_' + str(myRxn[i])
            LW = abs(float(myFlux[i])) * 0.005 # replacing 'str(5)' with 'LW' in the following line would yield a layover with the thickness proportional to the flux passing through the reactions. 
            myString += '\t' + reactionName + '\t' + str(5) + '\t' + color 
            rxnList.append(myString)
    
    
    OutputFileName = 'layover_%(someGene)s.txt' %{'someGene' : 'myModel'} #change outputFileName at convenience
    with open(OutputFileName, 'w') as file: 
        for k in range(len(rxnList)): 
            file.write(rxnList[k])
            file.write('\n')

    return 

if __name__ == '__main__': 
    # example of usage 
    import convertMyMedia
    import convertMyModel
    import cobra.io as co
    myModel = convertMyMedia.convertMyMedia(convertMyModel.convertMyModel('expression_Recon_ordered2.txt', co.read_sbml_model('Recon3.xml'))[0], 'minimal_media2.txt')
    makeLayover(myModel)
